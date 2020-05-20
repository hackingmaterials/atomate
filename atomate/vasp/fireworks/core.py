# coding: utf-8

import warnings
import copy

from atomate.vasp.config import (
    HALF_KPOINTS_FIRST_RELAX,
    RELAX_MAX_FORCE,
    VASP_CMD,
    DB_FILE,
    VDW_KERNEL_DIR,
)

"""
Defines standardized Fireworks that can be chained easily to perform various
sequences of VASP calculations.
"""

from fireworks import Firework

from pymatgen import Structure
from pymatgen.io.vasp.sets import (
    MPRelaxSet,
    MPScanRelaxSet,
    MITMDSet,
    MITRelaxSet,
    MPStaticSet,
    MPSOCSet,
)

from atomate.common.firetasks.glue_tasks import (
    PassCalcLocs,
    GzipDir,
    CopyFiles,
    DeleteFiles,
    CopyFilesFromCalcLoc,
)
from atomate.vasp.firetasks.glue_tasks import CopyVaspOutputs, pass_vasp_result
from atomate.vasp.firetasks.neb_tasks import TransferNEBTask
from atomate.vasp.firetasks.parse_outputs import VaspToDb, BoltztrapToDb
from atomate.vasp.firetasks.run_calc import (
    RunVaspCustodian,
    RunBoltztrap,
)
from atomate.vasp.firetasks.write_inputs import (
    WriteNormalmodeDisplacedPoscar,
    WriteTransmutedStructureIOSet,
    WriteVaspFromIOSet,
    WriteVaspHSEBSFromPrev,
    WriteVaspNSCFFromPrev,
    WriteVaspSOCFromPrev,
    WriteVaspStaticFromPrev,
    WriteVaspFromIOSetFromInterpolatedPOSCAR,
    UpdateScanRelaxBandgap,
    ModifyIncar,
)
from atomate.vasp.firetasks.neb_tasks import WriteNEBFromImages, WriteNEBFromEndpoints


class OptimizeFW(Firework):
    def __init__(
        self,
        structure,
        name="structure optimization",
        vasp_input_set=None,
        vasp_cmd=VASP_CMD,
        override_default_vasp_params=None,
        ediffg=None,
        db_file=DB_FILE,
        force_gamma=True,
        job_type="double_relaxation_run",
        max_force_threshold=RELAX_MAX_FORCE,
        auto_npar=">>auto_npar<<",
        half_kpts_first_relax=HALF_KPOINTS_FIRST_RELAX,
        parents=None,
        **kwargs
    ):
        """
        Optimize the given structure.

        Args:
            structure (Structure): Input structure.
            name (str): Name for the Firework.
            vasp_input_set (VaspInputSet): input set to use. Defaults to MPRelaxSet() if None.
            override_default_vasp_params (dict): If this is not None, these params are passed to
                the default vasp_input_set, i.e., MPRelaxSet. This allows one to easily override
                some settings, e.g., user_incar_settings, etc.
            vasp_cmd (str): Command to run vasp.
            ediffg (float): Shortcut to set ediffg in certain jobs
            db_file (str): Path to file specifying db credentials to place output parsing.
            force_gamma (bool): Force gamma centered kpoint generation
            job_type (str): custodian job type (default "double_relaxation_run")
            max_force_threshold (float): max force on a site allowed at end; otherwise, reject job
            auto_npar (bool or str): whether to set auto_npar. defaults to env_chk: ">>auto_npar<<"
            half_kpts_first_relax (bool): whether to use half the kpoints for the first relaxation
            parents ([Firework]): Parents of this particular Firework.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        override_default_vasp_params = override_default_vasp_params or {}
        vasp_input_set = vasp_input_set or MPRelaxSet(
            structure, force_gamma=force_gamma, **override_default_vasp_params
        )

        if vasp_input_set.incar["ISIF"] in (0, 1, 2, 7) and job_type == "double_relaxation":
            warnings.warn(
                "A double relaxation run might not be appropriate with ISIF {}".format(
                    vasp_input_set.incar["ISIF"]))
        
        t = []
        t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))
        t.append(
            RunVaspCustodian(
                vasp_cmd=vasp_cmd,
                job_type=job_type,
                max_force_threshold=max_force_threshold,
                ediffg=ediffg,
                auto_npar=auto_npar,
                half_kpts_first_relax=half_kpts_first_relax,
            )
        )
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, additional_fields={"task_label": name}))
        super(OptimizeFW, self).__init__(
            t,
            parents=parents,
            name="{}-{}".format(structure.composition.reduced_formula, name),
            **kwargs
        )


class ScanOptimizeFW(Firework):
    def __init__(
        self,
        structure,
        name="SCAN structure optimization",
        vasp_input_set=None,
        vasp_cmd=VASP_CMD,
        override_default_vasp_params=None,
        db_file=DB_FILE,
        vdw_kernel_dir=VDW_KERNEL_DIR,
        parents=None,
        **kwargs
    ):
        """
        Structure optimization using the SCAN metaGGA functional.

        This workflow performs a 3-step optmization. The first step ('relax1')
        is a conventional GGA run relaxation that initializes the geometry and
        calculates the bandgap of the structure. The bandgap is used to update 
        the KSPACING parameter, which sets the appropriate number of k-points 
        for the structure. The second step ('.relax2') is a static GGA 
        calculation that computes wavefunctions using the updated number of 
        k-points. The third step ('relax3') is a SCAN relaxation.

        By default, .relax1 and .relax2 are force converged with
        EDIFFG = -0.05, and .relax3 is force converged with EDIFFG=-0.02.

        Args:
            structure (Structure): Input structure.
            name (str): Name for the Firework.
            vasp_input_set (VaspInputSet): input set to use. Defaults to
                MPScanRelaxSet() if None.
            override_default_vasp_params (dict): If this is not None, and
                vasp_input_set is None, these params are passed to the default
                vasp_input_set, i.e., MPScanRelaxSet. This allows one to easily
                override some settings, e.g., bandgap, user_incar_settings, etc.
            vasp_cmd (str): Command to run vasp. Supports env_chk.
            vdw_kernel_dir (str): Directory containing the pre-compiled VdW
                kernel. Supports env_chk.
            db_file (str): Path to file specifying db credentials to place
                output parsing. Supports env_chk.
            parents ([Firework]): Parents of this particular Firework.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """
        override_default_vasp_params = override_default_vasp_params or {}
        orig_input_set = vasp_input_set or MPScanRelaxSet(
            structure, **override_default_vasp_params
        )

        # Raise a warning if the InputSet is not MPScanRelaxSet, because the
        # kspacing calculation from bandgap is only supported in MPScanRelaxSet.
        if not isinstance(orig_input_set, MPScanRelaxSet):
            raise UserWarning(
                "You have specified a vasp_input_set other than \
                               MPScanRelaxSet. Automatic adjustment of kspacing\
                               is not supported by this InputSet."
            )

        t = []
        # write the VASP input files based on MPScanRelaxSet
        t.append(WriteVaspFromIOSet(structure=structure, 
                                    vasp_input_set=orig_input_set
                                    )
                 )

        # pass the CalcLoc so that CopyFilesFromCalcLoc can find the directory
        t.append(PassCalcLocs(name=name))

        # Copy the pre-compiled VdW kernel for VASP, if required
        if orig_input_set.vdw is not None:
            t.append(CopyFiles(from_dir=vdw_kernel_dir))

        # Copy original inputs with the ".orig" suffix
        t.append(
            CopyFilesFromCalcLoc(
                calc_loc=True,
                name_append=".orig",
                exclude_files=["vdw_kernel.bindat", "FW.json", "FW--*"],
            )
        )

        # Update the INCAR for the GGA preconditioning step
        # Disable writing the WAVECAR because the no. of k-points will likely
        # change before the next step in the calculation
        pre_opt_settings = {"_set": {"METAGGA": None,
                                     "EDIFFG": -0.05,
                                     "LWAVE": False}}

        # Disable vdW for the precondition step
        if orig_input_set.incar.get("LUSE_VDW", None):
            pre_opt_settings.update({"_unset": {"LUSE_VDW": True,
                                                "BPARAM": 15.7}})

        t.append(ModifyIncar(incar_dictmod=pre_opt_settings))

        # Run the GGA .relax1 step
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, 
                                  job_type="normal_no_backup",
                                  gzip_output=False
                                  )
                 )

        # Copy GGA outputs with '.relax1' suffix
        # by the subsequent UpdateScanRelaxBandgap Firetask
        t.append(
            CopyFilesFromCalcLoc(
                calc_loc=True,
                name_append=".relax1",
                exclude_files=["vdw_kernel.bindat", "FW.json", "FW--*", "*.orig"],
            )
        )

        # Create a new InputSet and write new inputs based on the bandgap
        # set ICHARG and ISTART = 1 to start the calc from the previous charge
        # density and WAVECAR, respectively
        other_params = copy.deepcopy(orig_input_set.kwargs)
        if other_params.get("user_incar_settings"):
            other_params["user_incar_settings"]["ISTART"] = 1
            other_params["user_incar_settings"]["ICHARG"] = 1
        else:
            other_params["user_incar_settings"] = {"ISTART": 1,
                                                   "ICHARG": 1}
        t.append(UpdateScanRelaxBandgap(override_default_vasp_params=other_params))

        # Store the INCAR generated by UpdateScanRelaxBandGap for later use
        t.append(ModifyIncar(output_filename="INCAR.temp"))

        # Run a GGA static (.relax2) to initialize the wavefunction
        # In addition to the previous INCAR updates used in .relax1, output the
        # WAVECAR with LWAVE True and set ISTART to 0, since there is no WAVECAR
        # from the previous calculation
        pre_opt_settings2 = {"_set": {"METAGGA": None,
                                      "EDIFFG": -0.05,
                                      "LWAVE": True,
                                      "NSW": 0,
                                      "ISTART": 0}}
        if orig_input_set.incar.get("LUSE_VDW", None):
            pre_opt_settings2.update({"_unset": {"LUSE_VDW": True,
                                                 "BPARAM": 15.7}})

        # Update the INCAR for the GGA static run
        t.append(ModifyIncar(incar_dictmod=pre_opt_settings2))

        # Run the GGA static .relax2 step
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd,
                                  job_type="normal_no_backup",
                                  gzip_output=False
                                  )
                 )

        # Copy GGA outputs with '.relax2' suffix
        t.append(
            CopyFilesFromCalcLoc(
                calc_loc=True,
                name_append=".relax2",
                exclude_files=[
                    "vdw_kernel.bindat",
                    "FW.json",
                    "FW--*",
                    "*.orig",
                    "*.relax1",
                    "INCAR.temp"
                    ],
            )
        )

        # Reset the INCAR to the settings given by UpdateScanRelaxBandgap
        # by copying INCAR.temp to INCAR
        t.append(ModifyIncar(input_filename="INCAR.temp"))

        # Run the SCAN optimization step
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, job_type="normal_no_backup",
                                  gzip_output=False))

        # Copy SCAN outputs with '.relax3' suffix
        t.append(
            CopyFilesFromCalcLoc(
                calc_loc=True,
                name_append=".relax3",
                exclude_files=[
                    "vdw_kernel.bindat",
                    "FW.json",
                    "FW--*",
                    "*.orig",
                    "*.relax1",
                    "*.relax2",
                    "INCAR.temp"
                ],
            )
        )

        # Delete the VdW kernel, WAVECAR, custodian.json, and VASP output files
        # that have been copied to .relax3
        # Deleting custodian.json is necessary to avoid double-counting
        # custodian.json and custodian.json.relax3 when the output is parsed
        t.append(
            DeleteFiles(
                files=[
                    "vdw_kernel.bindat",
                    "WAVECAR*",  # All WAVECARs,
                    "custodian.json",
                    "*CAR",  # All files that end in "CAR"
                    "EIGENVAL",
                    "AECCAR?",
                    "IBZKPT",
                    "LOCPOT",
                    "REPORT",
                    "std_err.txt",
                    "vasp.out",
                    "CHG",
                    "PCDAT",
                    "vasprun.xml",
                    "INCAR.temp"
                ]
            )
        )

        # Parse the outputs into the database
        t.append(VaspToDb(db_file=db_file, additional_fields={"task_label": name}))

        # gzip the output
        t.append(GzipDir())

        super(ScanOptimizeFW, self).__init__(
            t,
            parents=parents,
            name="{}-{}".format(structure.composition.reduced_formula, name),
            **kwargs
        )


class StaticFW(Firework):
    def __init__(
        self,
        structure=None,
        name="static",
        vasp_input_set=None,
        vasp_input_set_params=None,
        vasp_cmd=VASP_CMD,
        prev_calc_loc=True,
        prev_calc_dir=None,
        db_file=DB_FILE,
        vasptodb_kwargs=None,
        parents=None,
        **kwargs
    ):
        """
        Standard static calculation Firework - either from a previous location or from a structure.

        Args:
            structure (Structure): Input structure. Note that for prev_calc_loc jobs, the structure
                is only used to set the name of the FW and any structure with the same composition
                can be used.
            name (str): Name for the Firework.
            vasp_input_set (VaspInputSet): input set to use (for jobs w/no parents)
                Defaults to MPStaticSet() if None.
            vasp_input_set_params (dict): Dict of vasp_input_set kwargs.
            vasp_cmd (str): Command to run vasp.
            prev_calc_loc (bool or str): If true (default), copies outputs from previous calc. If
                a str value, retrieves a previous calculation output by name. If False/None, will create
                new static calculation using the provided structure.
            prev_calc_dir (str): Path to a previous calculation to copy from
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            vasptodb_kwargs (dict): kwargs to pass to VaspToDb
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        t = []

        vasp_input_set_params = vasp_input_set_params or {}
        vasptodb_kwargs = vasptodb_kwargs or {}
        if "additional_fields" not in vasptodb_kwargs:
            vasptodb_kwargs["additional_fields"] = {}
        vasptodb_kwargs["additional_fields"]["task_label"] = name

        fw_name = "{}-{}".format(
            structure.composition.reduced_formula if structure else "unknown", name
        )

        if prev_calc_dir:
            t.append(CopyVaspOutputs(calc_dir=prev_calc_dir, contcar_to_poscar=True))
            t.append(WriteVaspStaticFromPrev(other_params=vasp_input_set_params))
        elif parents:
            if prev_calc_loc:
                t.append(
                    CopyVaspOutputs(calc_loc=prev_calc_loc, contcar_to_poscar=True)
                )
            t.append(WriteVaspStaticFromPrev(other_params=vasp_input_set_params))
        elif structure:
            vasp_input_set = vasp_input_set or MPStaticSet(
                structure, **vasp_input_set_params
            )
            t.append(
                WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set)
            )
        else:
            raise ValueError("Must specify structure or previous calculation")

        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, **vasptodb_kwargs))
        super(StaticFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)


class StaticInterpolateFW(Firework):
    def __init__(
        self,
        structure,
        start,
        end,
        name="static",
        vasp_input_set="MPStaticSet",
        vasp_input_set_params=None,
        vasp_cmd=VASP_CMD,
        db_file=DB_FILE,
        parents=None,
        this_image=None,
        nimages=None,
        autosort_tol=0,
        **kwargs
    ):
        """
        Standard static calculation Firework that interpolates structures from two previous calculations.

        Args:
            structure (Structure): Input structure used to name FireWork.
            start (str): PassCalcLoc name of StaticFW or RelaxFW run of starting structure.
            end (str): PassCalcLoc name of StaticFW or RelaxFW run of ending structure.
            name (str): Name for the Firework.
            vasp_input_set (str): Input set to use. Defaults to MPStaticSet.
            vasp_input_set_params (dict): Dict of vasp_input_set_kwargs.
            vasp_cmd (str): Command to run vasp.
            copy_vasp_outputs (bool): Whether to copy outputs from previous run. Defaults to True.
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            this_image (int): which interpolation to use for this run
            nimages (int): number of interpolations
            autosort_tol (float): a distance tolerance in angstrom in which
                to automatically sort end_structure to match to the closest
                points in this particular structure.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        t = []

        vasp_input_set_params = vasp_input_set_params or {}

        t.append(
            WriteVaspFromIOSetFromInterpolatedPOSCAR(
                start=start,
                end=end,
                this_image=this_image,
                nimages=nimages,
                autosort_tol=autosort_tol,
                vasp_input_set=vasp_input_set,
                vasp_input_params=vasp_input_set_params,
            )
        )

        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, additional_fields={"task_label": name}))

        super(StaticInterpolateFW, self).__init__(
            t,
            parents=parents,
            name="{}-{}".format(structure.composition.reduced_formula, name),
            **kwargs
        )


class HSEBSFW(Firework):
    def __init__(
        self,
        parents=None,
        prev_calc_dir=None,
        structure=None,
        mode="gap",
        name=None,
        vasp_cmd=VASP_CMD,
        db_file=DB_FILE,
        **kwargs
    ):
        """
        For getting a more accurate band gap or a full band structure with HSE - requires previous
        calculation that gives VBM/CBM info or the high-symmetry kpoints.

        Args:
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            prev_calc_dir (str): Path to a previous calculation to copy from
            structure (Structure): Input structure - used only to set the name of the FW.
            mode (string): options:
                "line" to get a full band structure along symmetry lines or
                "uniform" for uniform mesh band structure or
                "gap" to get the energy at the CBM and VBM
            name (str): Name for the Firework.
            vasp_cmd (str): Command to run vasp.
            db_file (str): Path to file specifying db credentials.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        name = name if name else "{} {}".format("hse", mode)

        fw_name = "{}-{}".format(
            structure.composition.reduced_formula if structure else "unknown", name
        )

        t = []
        if prev_calc_dir:
            t.append(
                CopyVaspOutputs(calc_dir=prev_calc_dir, additional_files=["CHGCAR"])
            )
        elif parents:
            t.append(CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"]))
        else:
            raise ValueError("Must specify a previous calculation for HSEBSFW")

        t.append(WriteVaspHSEBSFromPrev(prev_calc_dir=".", mode=mode))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd))
        t.append(PassCalcLocs(name=name))

        parse_dos = True if mode == "uniform" else False
        bandstructure_mode = mode if mode in ["line", "uniform"] else "line"

        t.append(
            VaspToDb(
                db_file=db_file,
                additional_fields={"task_label": name},
                parse_dos=parse_dos,
                bandstructure_mode=bandstructure_mode,
            )
        )
        super(HSEBSFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)


class NonSCFFW(Firework):
    def __init__(
        self,
        parents=None,
        prev_calc_dir=None,
        structure=None,
        name="nscf",
        mode="uniform",
        vasp_cmd=VASP_CMD,
        copy_vasp_outputs=True,
        db_file=DB_FILE,
        input_set_overrides=None,
        **kwargs
    ):
        """
        Standard NonSCF Calculation Firework supporting uniform and line modes.

        Args:
            structure (Structure): Input structure - used only to set the name
                of the FW.
            name (str): Name for the Firework.
            mode (str): "uniform" or "line" mode.
            vasp_cmd (str): Command to run vasp.
            copy_vasp_outputs (bool): Whether to copy outputs from previous
                run. Defaults to True.
            prev_calc_dir (str): Path to a previous calculation to copy from
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework.
                FW or list of FWS.
            input_set_overrides (dict): Arguments passed to the
                "from_prev_calc" method of the MPNonSCFSet. This parameter
                allows a user to modify the default values of the input set.
                For example, passing the key value pair
                    {'reciprocal_density': 1000}
                will override default k-point meshes for uniform calculations.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        input_set_overrides = input_set_overrides or {}

        fw_name = "{}-{} {}".format(
            structure.composition.reduced_formula if structure else "unknown",
            name,
            mode,
        )
        t = []

        if prev_calc_dir:
            t.append(
                CopyVaspOutputs(calc_dir=prev_calc_dir, additional_files=["CHGCAR"])
            )
        elif parents:
            t.append(CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"]))
        else:
            raise ValueError("Must specify previous calculation for NonSCFFW")

        mode = mode.lower()
        if mode == "uniform":
            t.append(
                WriteVaspNSCFFromPrev(
                    prev_calc_dir=".", mode="uniform", **input_set_overrides
                )
            )
        else:
            t.append(
                WriteVaspNSCFFromPrev(
                    prev_calc_dir=".", mode="line", **input_set_overrides
                )
            )

        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"))
        t.append(PassCalcLocs(name=name))
        t.append(
            VaspToDb(
                db_file=db_file,
                additional_fields={"task_label": name + " " + mode},
                parse_dos=(mode == "uniform"),
                bandstructure_mode=mode,
            )
        )

        super(NonSCFFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)


class LepsFW(Firework):
    def __init__(
        self,
        structure,
        name="static dielectric",
        vasp_cmd=VASP_CMD,
        copy_vasp_outputs=True,
        db_file=DB_FILE,
        parents=None,
        phonon=False,
        mode=None,
        displacement=None,
        user_incar_settings=None,
        **kwargs
    ):
        """
        Standard static calculation Firework for dielectric constants using DFPT.

        Args:
            structure (Structure): Input structure. If copy_vasp_outputs, used only to set the
                name of the FW.
            name (str): Name for the Firework.
            vasp_cmd (str): Command to run vasp.
            copy_vasp_outputs (bool): Whether to copy outputs from previous
                run. Defaults to True.
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework.
                FW or list of FWS.
            phonon (bool): Whether or not to extract normal modes and pass it. This argument along
                with the mode and displacement arguments must be set for the calculation of
                dielectric constant in the Raman tensor workflow.
            mode (int): normal mode index.
            displacement (float): displacement along the normal mode in Angstroms.
            user_incar_settings (dict): Parameters in INCAR to override
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        warnings.warn(
            "This firework will be removed soon. Use DFPTFW and/or RamanFW fireworks."
        )
        user_incar_settings = user_incar_settings or {}
        t = []

        if copy_vasp_outputs:
            t.append(
                CopyVaspOutputs(
                    calc_loc=True, additional_files=["CHGCAR"], contcar_to_poscar=True
                )
            )
            t.append(
                WriteVaspStaticFromPrev(
                    lepsilon=True,
                    other_params={"user_incar_settings": user_incar_settings},
                )
            )
        else:
            vasp_input_set = MPStaticSet(
                structure, lepsilon=True, user_incar_settings=user_incar_settings
            )
            t.append(
                WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set)
            )

        if phonon:
            if mode is None and displacement is None:
                name = "{} {}".format("phonon", name)
                t.append(RunVaspCustodian(vasp_cmd=vasp_cmd))
                t.append(
                    pass_vasp_result(
                        {
                            "structure": "a>>final_structure",
                            "eigenvals": "a>>normalmode_eigenvals",
                            "eigenvecs": "a>>normalmode_eigenvecs",
                        },
                        parse_eigen=True,
                        mod_spec_key="normalmodes",
                    )
                )
            else:
                name = "raman_{}_{} {}".format(str(mode), str(displacement), name)
                key = (
                    "{}_{}".format(mode, displacement)
                    .replace("-", "m")
                    .replace(".", "d")
                )
                pass_fw = pass_vasp_result(
                    pass_dict={
                        "mode": mode,
                        "displacement": displacement,
                        "epsilon": "a>>epsilon_static",
                    },
                    mod_spec_key="raman_epsilon->" + key,
                    parse_eigen=True,
                )
                t.extend(
                    [
                        WriteNormalmodeDisplacedPoscar(
                            mode=mode, displacement=displacement
                        ),
                        RunVaspCustodian(vasp_cmd=vasp_cmd),
                        pass_fw,
                    ]
                )
        else:
            t.append(RunVaspCustodian(vasp_cmd=vasp_cmd))

        t.extend(
            [
                PassCalcLocs(name=name),
                VaspToDb(db_file=db_file, additional_fields={"task_label": name}),
            ]
        )

        super(LepsFW, self).__init__(
            t,
            parents=parents,
            name="{}-{}".format(structure.composition.reduced_formula, name),
            **kwargs
        )


class DFPTFW(Firework):
    def __init__(
        self,
        structure=None,
        prev_calc_dir=None,
        name="static dielectric",
        vasp_cmd=VASP_CMD,
        copy_vasp_outputs=True,
        lepsilon=True,
        db_file=DB_FILE,
        parents=None,
        user_incar_settings=None,
        pass_nm_results=False,
        **kwargs
    ):
        """
         Static DFPT calculation Firework

        Args:
            structure (Structure): Input structure. If copy_vasp_outputs, used only to set the
                name of the FW.
            name (str): Name for the Firework.
            lepsilon (bool): Turn on LEPSILON to calculate polar properties
            vasp_cmd (str): Command to run vasp.
            copy_vasp_outputs (str or bool): Whether to copy outputs from previous
                run. Defaults to True.
            prev_calc_dir (str): Path to a previous calculation to copy from
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework.
                FW or list of FWS.
            user_incar_settings (dict): Parameters in INCAR to override
            pass_nm_results (bool): if true the normal mode eigen vals and vecs are passed so that
                next firework can use it.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        name = "static dielectric" if lepsilon else "phonon"

        fw_name = "{}-{}".format(
            structure.composition.reduced_formula if structure else "unknown", name
        )

        user_incar_settings = user_incar_settings or {}
        t = []

        if prev_calc_dir:
            t.append(CopyVaspOutputs(calc_dir=prev_calc_dir, contcar_to_poscar=True))
            t.append(
                WriteVaspStaticFromPrev(
                    lepsilon=lepsilon,
                    other_params={
                        "user_incar_settings": user_incar_settings,
                        "force_gamma": True,
                    },
                )
            )
        elif parents and copy_vasp_outputs:
            t.append(CopyVaspOutputs(calc_loc=True, contcar_to_poscar=True))
            t.append(
                WriteVaspStaticFromPrev(
                    lepsilon=lepsilon,
                    other_params={
                        "user_incar_settings": user_incar_settings,
                        "force_gamma": True,
                    },
                )
            )
        elif structure:
            vasp_input_set = MPStaticSet(
                structure,
                lepsilon=lepsilon,
                force_gamma=True,
                user_incar_settings=user_incar_settings,
            )
            t.append(
                WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set)
            )
        else:
            raise ValueError("Must specify structure or previous calculation")

        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd))

        if pass_nm_results:
            t.append(
                pass_vasp_result(
                    {
                        "structure": "a>>final_structure",
                        "eigenvals": "a>>normalmode_eigenvals",
                        "eigenvecs": "a>>normalmode_eigenvecs",
                    },
                    parse_eigen=True,
                    mod_spec_key="normalmodes",
                )
            )

        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, additional_fields={"task_label": name}))

        super(DFPTFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)


class RamanFW(Firework):
    def __init__(
        self,
        mode,
        displacement,
        prev_calc_dir=None,
        structure=None,
        name="raman",
        vasp_cmd=VASP_CMD,
        db_file=DB_FILE,
        parents=None,
        user_incar_settings=None,
        **kwargs
    ):
        """
        Static calculation Firework that computes the DFPT dielectric constant for
        structure displaced along the given normal mode direction.

        Args:
            structure (Structure): Input structure. If copy_vasp_outputs, used only to set the
                name of the FW.
            mode (int): normal mode index.
            displacement (float): displacement along the normal mode in Angstroms.
            name (str): Name for the Firework.
            prev_calc_dir (str): Path to a previous calculation to copy from
            vasp_cmd (str): Command to run vasp.
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework.
                FW or list of FWS.
            user_incar_settings (dict): Parameters in INCAR to override
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        name = "{}_{}_{}".format(name, str(mode), str(displacement))
        fw_name = "{}-{}".format(
            structure.composition.reduced_formula if structure else "unknown", name
        )

        user_incar_settings = user_incar_settings or {}

        t = []

        if prev_calc_dir:
            t.append(CopyVaspOutputs(calc_dir=prev_calc_dir, contcar_to_poscar=True))
        elif parents:
            t.append(CopyVaspOutputs(calc_loc=True, contcar_to_poscar=True))
        else:
            raise ValueError("Must specify a previous calculation")

        t.append(
            WriteVaspStaticFromPrev(
                lepsilon=True, other_params={"user_incar_settings": user_incar_settings}
            )
        )

        t.append(WriteNormalmodeDisplacedPoscar(mode=mode, displacement=displacement))

        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd))

        key = "{}_{}".format(mode, displacement).replace("-", "m").replace(".", "d")
        t.append(
            pass_vasp_result(
                pass_dict={
                    "mode": mode,
                    "displacement": displacement,
                    "epsilon": "a>>epsilon_static",
                },
                mod_spec_key="raman_epsilon->{}".format(key),
                parse_eigen=True,
            )
        )

        t.append(PassCalcLocs(name=name))

        t.append(VaspToDb(db_file=db_file, additional_fields={"task_label": name}))

        super(RamanFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)


class SOCFW(Firework):
    def __init__(
        self,
        magmom,
        structure=None,
        name="spin-orbit coupling",
        saxis=(0, 0, 1),
        prev_calc_dir=None,
        vasp_cmd="vasp_ncl",
        copy_vasp_outputs=True,
        db_file=None,
        parents=None,
        **kwargs
    ):
        """
        Firework for spin orbit coupling calculation.

        Args:
            structure (Structure): Input structure. If copy_vasp_outputs, used only to set the
                name of the FW.
            name (str): Name for the Firework.
            prev_calc_dir (str): Path to a previous calculation to copy from
            vasp_cmd (str): Command to run vasp.
            copy_vasp_outputs (bool): Whether to copy outputs from previous
                run. Defaults to True.
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework.
                FW or list of FWS.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        fw_name = "{}-{}".format(
            structure.composition.reduced_formula if structure else "unknown", name
        )

        t = []
        if prev_calc_dir:
            t.append(
                CopyVaspOutputs(
                    calc_dir=prev_calc_dir,
                    additional_files=["CHGCAR"],
                    contcar_to_poscar=True,
                )
            )
            t.append(
                WriteVaspSOCFromPrev(prev_calc_dir=".", magmom=magmom, saxis=saxis)
            )
        elif parents and copy_vasp_outputs:
            t.append(
                CopyVaspOutputs(
                    calc_loc=True, additional_files=["CHGCAR"], contcar_to_poscar=True
                )
            )
            t.append(
                WriteVaspSOCFromPrev(prev_calc_dir=".", magmom=magmom, saxis=saxis)
            )
        elif structure:
            vasp_input_set = MPSOCSet(structure)
            t.append(
                WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set)
            )
        else:
            raise ValueError("Must specify structure or previous calculation.")

        t.extend(
            [
                RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"),
                PassCalcLocs(name=name),
                VaspToDb(db_file=db_file, additional_fields={"task_label": name}),
            ]
        )
        super(SOCFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)


class TransmuterFW(Firework):
    def __init__(
        self,
        structure,
        transformations,
        transformation_params=None,
        vasp_input_set=None,
        prev_calc_dir=None,
        name="structure transmuter",
        vasp_cmd=VASP_CMD,
        copy_vasp_outputs=True,
        db_file=DB_FILE,
        parents=None,
        override_default_vasp_params=None,
        **kwargs
    ):
        """
        Apply the transformations to the input structure, write the input set corresponding
        to the transformed structure, and run vasp on them.  Note that if a transformation yields
        many structures from one, only the last structure in the list is used.

        Args:
            structure (Structure): Input structure.
            transformations (list): list of names of transformation classes as defined in
                the modules in pymatgen.transformations.
                eg:  transformations=['DeformStructureTransformation', 'SupercellTransformation']
            transformation_params (list): list of dicts where each dict specify the input
                parameters to instantiate the transformation class in the transformations list.
            vasp_input_set (VaspInputSet): VASP input set, used to write the input set for the
                transmuted structure.
            name (string): Name for the Firework.
            vasp_cmd (string): Command to run vasp.
            copy_vasp_outputs (bool): Whether to copy outputs from previous run. Defaults to True.
            prev_calc_dir (str): Path to a previous calculation to copy from
            db_file (string): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            override_default_vasp_params (dict): additional user input settings for vasp_input_set.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        fw_name = "{}-{}".format(structure.composition.reduced_formula, name)
        override_default_vasp_params = override_default_vasp_params or {}
        t = []

        vasp_input_set = vasp_input_set or MPStaticSet(
            structure, force_gamma=True, **override_default_vasp_params
        )

        if prev_calc_dir:
            t.append(CopyVaspOutputs(calc_dir=prev_calc_dir, contcar_to_poscar=True))
            t.append(
                WriteTransmutedStructureIOSet(
                    transformations=transformations,
                    transformation_params=transformation_params,
                    vasp_input_set=vasp_input_set,
                    override_default_vasp_params=override_default_vasp_params,
                    prev_calc_dir=".",
                )
            )
        elif copy_vasp_outputs:
            t.append(CopyVaspOutputs(calc_loc=True, contcar_to_poscar=True))
            t.append(
                WriteTransmutedStructureIOSet(
                    structure=structure,
                    transformations=transformations,
                    transformation_params=transformation_params,
                    vasp_input_set=vasp_input_set,
                    override_default_vasp_params=override_default_vasp_params,
                    prev_calc_dir=".",
                )
            )
        elif structure:
            t.append(
                WriteTransmutedStructureIOSet(
                    structure=structure,
                    transformations=transformations,
                    transformation_params=transformation_params,
                    vasp_input_set=vasp_input_set,
                    override_default_vasp_params=override_default_vasp_params,
                )
            )
        else:
            raise ValueError("Must specify structure or previous calculation")

        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd))
        t.append(PassCalcLocs(name=name))
        t.append(
            VaspToDb(
                db_file=db_file,
                additional_fields={
                    "task_label": name,
                    "transmuter": {
                        "transformations": transformations,
                        "transformation_params": transformation_params,
                    },
                },
            )
        )

        super(TransmuterFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)


class MDFW(Firework):
    def __init__(
        self,
        structure,
        start_temp,
        end_temp,
        nsteps,
        name="molecular dynamics",
        vasp_input_set=None,
        vasp_cmd=VASP_CMD,
        override_default_vasp_params=None,
        wall_time=19200,
        db_file=DB_FILE,
        parents=None,
        copy_vasp_outputs=True,
        **kwargs
    ):
        """
        Standard firework for a single MD run.

        Args:
            structure (Structure): Input structure.
            start_temp (float): Start temperature of MD run.
            end_temp (float): End temperature of MD run.
            nsteps (int): Number of MD steps
            name (string): Name for the Firework.
            vasp_input_set (string): string name for the VASP input set (e.g.,
                "MITMDVaspInputSet").
            vasp_cmd (string): Command to run vasp.
            override_default_vasp_params (dict): If this is not None,
                these params are passed to the default vasp_input_set, i.e.,
                MITMDSet. This allows one to easily override some
                settings, e.g., user_incar_settings, etc. Particular to MD,
                one can control time_step and all other settings of the input set.
            wall_time (int): Total wall time in seconds before writing STOPCAR.
            copy_vasp_outputs (bool): Whether to copy outputs from previous run. Defaults to True.
            db_file (string): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        override_default_vasp_params = override_default_vasp_params or {}
        vasp_input_set = vasp_input_set or MITMDSet(
            structure,
            start_temp=start_temp,
            end_temp=end_temp,
            nsteps=nsteps,
            **override_default_vasp_params
        )

        t = []
        if copy_vasp_outputs:
            t.append(
                CopyVaspOutputs(
                    calc_loc=True, additional_files=["CHGCAR"], contcar_to_poscar=True
                )
            )

        t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))
        t.append(
            RunVaspCustodian(
                vasp_cmd=vasp_cmd,
                gamma_vasp_cmd=">>gamma_vasp_cmd<<",
                handler_group="md",
                wall_time=wall_time,
            )
        )
        t.append(PassCalcLocs(name=name))
        t.append(
            VaspToDb(
                db_file=db_file,
                additional_fields={"task_label": name},
                defuse_unsuccessful=False,
            )
        )
        super(MDFW, self).__init__(
            t,
            parents=parents,
            name="{}-{}".format(structure.composition.reduced_formula, name),
            **kwargs
        )


class BoltztrapFW(Firework):
    def __init__(
        self,
        parents=None,
        structure=None,
        name="boltztrap",
        db_file=None,
        scissor=0.0,
        doping=None,
        tmax=1300,
        tgrid=50,
        prev_calc_dir=None,
        soc=False,
        additional_fields=None,
        **kwargs
    ):
        """
        Run Boltztrap (which includes writing bolztrap input files and parsing outputs). Assumes
        you have a previous FW with the calc_locs passed into the current FW.

        Args:
            structure (Structure): - only used for setting name of FW
            name (str): name of this FW
            db_file (str): path to the db file
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            prev_calc_dir (str): Path to a previous calculation to copy from
            scissor (float): if scissor > 0, apply scissor on the band structure so that new
                band gap = scissor (in eV)
            doping: ([float]) doping levels you want to compute
            tmax: (float) max temperature to evaluate
            tgrid: (float) temperature interval
            soc (bool): whether the band structure is calculated with spin-orbit coupling
            additional_fields (dict): fields added to the document such as user-defined tags or name, ids, etc
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        fw_name = "{}-{}".format(
            structure.composition.reduced_formula if structure else "unknown", name
        )

        additional_fields = additional_fields or {}

        t = []

        if prev_calc_dir:
            t.append(CopyVaspOutputs(calc_dir=prev_calc_dir, contcar_to_poscar=True))
        elif parents:
            t.append(CopyVaspOutputs(calc_loc=True, contcar_to_poscar=True))
        else:
            raise ValueError("Must specify structure or previous calculation")

        t.extend(
            [
                RunBoltztrap(
                    scissor=scissor, soc=soc, doping=doping, tmax=tmax, tgrid=tgrid
                ),
                BoltztrapToDb(db_file=db_file, additional_fields=additional_fields),
                PassCalcLocs(name=name),
            ]
        )

        super(BoltztrapFW, self).__init__(t, parents=parents, name=fw_name, **kwargs)


class NEBRelaxationFW(Firework):
    """
    Relaxation Firework in NEB Workflow.

    Task 1) Read in a structure with "st_label" ("rlx", "ep0" or "ep1") and generates input sets.
    Task 2) Run VASP using Custodian
    Task 3) Update structure to spec
    Task 4) Pass CalcLocs named "{}_dir".format(st_label)
    """

    def __init__(
        self,
        spec,
        label,
        user_incar_settings=None,
        user_kpoints_settings=None,
        additional_cust_args=None,
        **kwargs
    ):
        """
        Args:
            spec (dict): Specification of the job to run.
            label (str): "parent", "ep0" or "ep1"
            vasp_input_set (VaspInputSet): Input set to use.
            user_kpoints_settings (dict): Additional KPOINTS settings.
            additional_cust_args (dict): Other kwargs that are passed to RunVaspCustodian.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """

        # Get structure from spec
        assert label in ["parent", "ep0", "ep1"]
        structure_dict = spec[label]
        structure = Structure.from_dict(structure_dict)

        user_incar_settings = user_incar_settings or {}
        user_kpoints_settings = user_kpoints_settings or {}
        additional_cust_args = additional_cust_args or {}

        # Task 1: Write input sets
        if label == "parent":
            vasp_input_set = MITRelaxSet(
                structure,
                user_incar_settings=user_incar_settings,
                user_kpoints_settings=user_kpoints_settings,
            )
        else:  # label == "ep0" or "ep1"
            from pymatgen_diffusion.neb.io import MVLCINEBEndPointSet

            vasp_input_set = MVLCINEBEndPointSet(
                structure,
                user_incar_settings=user_incar_settings,
                user_kpoints_settings=user_kpoints_settings,
            )

        write_ep_task = WriteVaspFromIOSet(
            structure=structure, vasp_input_set=vasp_input_set
        )

        # Task 2: Run VASP using Custodian
        cust_args = {
            "job_type": "normal",
            "gzip_output": False,
            "handler_group": "no_handler",
        }
        cust_args.update(additional_cust_args)
        run_vasp = RunVaspCustodian(
            vasp_cmd=">>vasp_cmd<<", gamma_vasp_cmd=">>gamma_vasp_cmd<<", **cust_args
        )

        # Task 3, 4: Transfer and PassCalLocs
        tasks = [
            write_ep_task,
            run_vasp,
            TransferNEBTask(label=label),
            PassCalcLocs(name=label),
        ]

        super(NEBRelaxationFW, self).__init__(tasks, spec=spec, name=label, **kwargs)


class NEBFW(Firework):
    """
    CI-NEB Firework in NEB Workflow.

    Task 1) Read in image structures from spec and generates input sets.
            The group of structures are labeled with neb_label (1, 2...)
    Task 2) Run NEB VASP using Custodian
    Task 3) Update structure to spec
    Task 4) Pass CalcLocs named "neb_{}".format(neb_label)
    """

    def __init__(
        self,
        spec,
        neb_label,
        from_images=True,
        user_incar_settings=None,
        user_kpoints_settings=None,
        additional_cust_args=None,
        **kwargs
    ):
        """
        Args:
            spec (dict): Specification of the job to run.
            neb_label (str): "1", "2"..., label neb run.
            from_images (bool): Set True to initialize from image structures, False starting
                        from relaxed endpoint structures.
            user_incar_settings (dict): Additional INCAR settings.
            user_kpoints_settings (dict): Additional KPOINTS settings.
            additional_cust_args (dict): Other kwargs that are passed to RunVaspCustodian.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        assert neb_label.isdigit() and int(neb_label) >= 1
        label = "neb{}".format(neb_label)
        sort_tol = spec["sort_tol"]
        d_img = spec["d_img"]
        interpolation_type = spec["interpolation_type"]

        # Task 1: Write NEB input sets
        user_incar_settings = user_incar_settings or {}
        user_kpoints_settings = user_kpoints_settings or {}
        additional_cust_args = additional_cust_args or {}

        if from_images:
            write_neb_task = WriteNEBFromImages(
                neb_label=neb_label,
                user_incar_settings=user_incar_settings,
                user_kpoints_settings=user_kpoints_settings,
            )

        else:  # from endpoints
            write_neb_task = WriteNEBFromEndpoints(
                user_incar_settings=user_incar_settings,
                user_kpoints_settings=user_kpoints_settings,
                sort_tol=sort_tol,
                d_img=d_img,
                interpolation_type=interpolation_type,
            )

        # Task 2: Run NEB using Custodian
        cust_args = {
            "job_type": "neb",
            "gzip_output": False,
            "handler_group": "no_handler",
        }
        cust_args.update(additional_cust_args)
        run_neb_task = RunVaspCustodian(
            vasp_cmd=">>vasp_cmd<<", gamma_vasp_cmd=">>gamma_vasp_cmd<<", **cust_args
        )

        # Task 3, 4: Transfer and PassCalcLocs
        tasks = [
            write_neb_task,
            run_neb_task,
            TransferNEBTask(label=label),
            PassCalcLocs(name=label),
        ]

        super(NEBFW, self).__init__(tasks, spec=spec, name=label, **kwargs)
