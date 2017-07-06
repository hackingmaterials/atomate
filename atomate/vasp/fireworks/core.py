# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
Defines standardized Fireworks that can be chained easily to perform various
sequences of VASP calculations.
"""

from fireworks import Firework

from pymatgen import Structure
from pymatgen.io.vasp.sets import MPRelaxSet, MITMDSet, MITRelaxSet, MPStaticSet, MPSOCSet

from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.vasp.firetasks.glue_tasks import CopyVaspOutputs, pass_vasp_result
from atomate.vasp.firetasks.neb_tasks import TransferNEBTask
from atomate.vasp.firetasks.parse_outputs import VaspToDb, BoltztrapToDb
from atomate.vasp.firetasks.run_calc import RunVaspCustodian, RunBoltztrap
from atomate.vasp.firetasks.write_inputs import WriteNormalmodeDisplacedPoscar, \
    WriteTransmutedStructureIOSet, WriteVaspFromIOSet, WriteVaspHSEBSFromPrev, \
    WriteVaspNSCFFromPrev, WriteVaspSOCFromPrev, WriteVaspStaticFromPrev
from atomate.vasp.firetasks.neb_tasks import WriteNEBFromImages, WriteNEBFromEndpoints


class OptimizeFW(Firework):
    def __init__(self, structure, name="structure optimization", vasp_input_set=None,
                 vasp_cmd="vasp", override_default_vasp_params=None, ediffg=None, db_file=None,
                 force_gamma=True, job_type="double_relaxation_run", max_force_threshold=0.25,
                 auto_npar=">>auto_npar<<", parents=None, **kwargs):
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
            parents ([Firework]): Parents of this particular Firework.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        override_default_vasp_params = override_default_vasp_params or {}
        vasp_input_set = vasp_input_set or MPRelaxSet(structure, force_gamma=force_gamma,
                                                      **override_default_vasp_params)

        t = []
        t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, job_type=job_type,
                                  max_force_threshold=max_force_threshold, ediffg=ediffg,
                                  auto_npar=auto_npar))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, additional_fields={"task_label": name}))
        super(OptimizeFW, self).__init__(t, parents=parents, name="{}-{}".
                                         format(structure.composition.reduced_formula, name),
                                         **kwargs)


class StaticFW(Firework):
    def __init__(self, structure, name="static", vasp_input_set=None, vasp_cmd="vasp",
                 prev_calc_loc=True, db_file=None, parents=None, **kwargs):
        """
        Standard static calculation Firework - either from a previous location or from a structure.

        Args:
            structure (Structure): Input structure. Note that for prev_calc_loc jobs, the structure 
                is only used to set the name of the FW and any structure with the same composition 
                can be used.
            name (str): Name for the Firework.
            vasp_input_set (VaspInputSet): input set to use (for jobs w/no parents)
                Defaults to MPStaticSet() if None.
            vasp_cmd (str): Command to run vasp.
            prev_calc_loc (bool or str): If true (default), copies outputs from previous calc. If 
                a str value, grabs a previous calculation output by name. If False/None, will create
                new static calculation using the provided structure.
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """

        # TODO: @computron - I really don't like how you need to set the structure even for
        # prev_calc_loc jobs. Sometimes it makes appending new FWs to an existing workflow
        # difficult. Maybe think about how to remove this need? -computron

        t = []

        if parents:
            if prev_calc_loc:
                t.append(CopyVaspOutputs(calc_loc=prev_calc_loc, contcar_to_poscar=True))
            t.append(WriteVaspStaticFromPrev())
        else:
            vasp_input_set = vasp_input_set or MPStaticSet(structure)
            t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))

        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, additional_fields={"task_label": name}))
        super(StaticFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)


class HSEBSFW(Firework):
    def __init__(self, structure, parents, mode="gap", name=None, vasp_cmd="vasp", db_file=None,
                 **kwargs):
        """
        For getting a more accurate band gap or a full band structure with HSE - requires previous
        calculation that gives VBM/CBM info or the high-symmetry kpoints.

        Args:
            structure (Structure): Input structure - used only to set the name of the FW.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
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

        t = []
        t.append(CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"]))
        t.append(WriteVaspHSEBSFromPrev(prev_calc_dir='.', mode=mode))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, additional_fields={"task_label": name}))
        super(HSEBSFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)


class NonSCFFW(Firework):
    def __init__(self, structure, name="nscf", mode="uniform", vasp_cmd="vasp",
                 copy_vasp_outputs=True, db_file=None, parents=None, **kwargs):
        """
        Standard NonSCF Calculation Firework supporting both
        uniform and line modes.

        Args:
            structure (Structure): Input structure - used only to set the name of the FW.
            name (str): Name for the Firework.
            mode (str): "uniform" or "line" mode.
            vasp_cmd (str): Command to run vasp.
            copy_vasp_outputs (bool): Whether to copy outputs from previous
                run. Defaults to True.
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework.
                FW or list of FWS.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        t = []
        if copy_vasp_outputs:
            t.append(CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"]))
        mode = mode.lower()
        if mode == "uniform":
            t.append(
                WriteVaspNSCFFromPrev(prev_calc_dir=".", mode="uniform", reciprocal_density=1000))
        else:
            t.append(WriteVaspNSCFFromPrev(prev_calc_dir=".", mode="line", reciprocal_density=20))

        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, additional_fields={"task_label": name + " " + mode},
                          parse_dos=(mode == "uniform"), bandstructure_mode=mode))

        super(NonSCFFW, self).__init__(t, parents=parents, name="%s-%s %s" % (
            structure.composition.reduced_formula, name, mode), **kwargs)


class LepsFW(Firework):
    def __init__(self, structure, name="static dielectric", vasp_cmd="vasp", copy_vasp_outputs=True,
                 db_file=None, parents=None, phonon=False, mode=None, displacement=None,
                 user_incar_settings=None, **kwargs):
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
        user_incar_settings = user_incar_settings or {}
        t = []

        if copy_vasp_outputs:
            t.append(CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"],
                                     contcar_to_poscar=True))
            t.append(WriteVaspStaticFromPrev(lepsilon=True,
                                             other_params={
                                                 'user_incar_settings': user_incar_settings}))
        else:
            vasp_input_set = MPStaticSet(structure, lepsilon=True,
                                         user_incar_settings=user_incar_settings)
            t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))

        if phonon:
            if mode is None and displacement is None:
                # TODO: @matk86 - I really don't understand the point of setting phonon=True w/o
                # the mode or displacement. It seems to simply be running a regular static run and
                # doing nothing else than changing the name of the Firework and perhaps passing
                # normal modes data (but it's unclear how that data is generated. Why would anyone
                # want to do this? -computron
                name = "{} {}".format("phonon", name)

                # TODO: @matk86 - not sure why this line is here. I understand you are trying to
                # keep the code short but the logic is very confusing. The RunVaspCustodian is
                # common to all 3 situations (phonon=F, phonon=T/mode=F, phonon=T/mode=T) yet is
                # duplicated in each place. Instead, for better clarity construct the Firework
                # sequentially (write inputs, run vasp, parse output data, pass output data) and
                # add if/else for phonons where it is needed. Any name overrides can go near the
                # top of the Firework. -computron
                t.append(RunVaspCustodian(vasp_cmd=vasp_cmd))
                t.append(pass_vasp_result({"structure": "a>>final_structure",
                                           "eigenvals": "a>>normalmode_eigenvals",
                                           "eigenvecs": "a>>normalmode_eigenvecs"}, parse_eigen=True,
                                           mod_spec_key="normalmodes"))
            else:
                # TODO: @matk86 - Why is this second calculation being tacked on to the first one?
                # It looks like it will overwrite INCAR/CHGCAR/etc of the first calculation and
                # thus completely remove access to the original output files. Shouldn't this be a
                # new Firework rather than a second calculation in the same Firework? -computron
                name = "raman_{}_{} {}".format(str(mode), str(displacement), name)
                key = "{}_{}".format(mode, displacement).replace('-', 'm').replace('.', 'd')
                pass_fw = pass_vasp_result(pass_dict={"mode": mode, "displacement": displacement,
                                                      "epsilon": "a>>epsilon_static"},
                                           mod_spec_key="raman_epsilon->" + key,
                                           parse_eigen=True)
                t.extend([WriteNormalmodeDisplacedPoscar(mode=mode, displacement=displacement),
                          RunVaspCustodian(vasp_cmd=vasp_cmd), pass_fw])
        else:
            t.append(RunVaspCustodian(vasp_cmd=vasp_cmd))

        t.extend([PassCalcLocs(name=name),
                  VaspToDb(db_file=db_file, additional_fields={"task_label": name})])

        super(LepsFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)


class SOCFW(Firework):
    def __init__(self, structure, magmom, name="spin-orbit coupling", saxis=(0, 0, 1),
                 vasp_cmd="vasp_ncl", copy_vasp_outputs=True, db_file=None, parents=None, **kwargs):
        """
        Firework for spin orbit coupling calculation.

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
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        t = []

        if copy_vasp_outputs:
            t.append(CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"],
                                     contcar_to_poscar=True))
            t.append(WriteVaspSOCFromPrev(prev_calc_dir=".", magmom=magmom, saxis=saxis))
        else:
            vasp_input_set = MPSOCSet(structure)
            t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))

        t.extend([RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"),
                  PassCalcLocs(name=name),
                  VaspToDb(db_file=db_file, additional_fields={"task_label": name})])
        super(SOCFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)


class TransmuterFW(Firework):
    def __init__(self, structure, transformations, transformation_params=None, vasp_input_set=None,
                 name="structure transmuter", vasp_cmd="vasp", copy_vasp_outputs=True, db_file=None,
                 parents=None, override_default_vasp_params={}, **kwargs):
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
            db_file (string): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            override_default_vasp_params (dict): additional user input settings for vasp_input_set.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        t = []

        vasp_input_set = vasp_input_set or MPStaticSet(structure, force_gamma=True,
                                                       **override_default_vasp_params)

        if copy_vasp_outputs:
            t.append(CopyVaspOutputs(calc_loc=True, contcar_to_poscar=True))

        prev_calc_dir = "." if copy_vasp_outputs else None
        t.append(
            WriteTransmutedStructureIOSet(structure=structure, transformations=transformations,
                                          transformation_params=transformation_params,
                                          vasp_input_set=vasp_input_set,
                                          override_default_vasp_params=override_default_vasp_params,
                                          prev_calc_dir=prev_calc_dir))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file,
                          additional_fields={
                                  "task_label": name,
                                  "transmuter": {"transformations": transformations,
                                                 "transformation_params": transformation_params}
                              }))

        super(TransmuterFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)


class MDFW(Firework):
    def __init__(self, structure, start_temp, end_temp, nsteps, name="molecular dynamics",
                 vasp_input_set=None, vasp_cmd="vasp", override_default_vasp_params=None,
                 wall_time=19200, db_file=None, parents=None, copy_vasp_outputs=True, **kwargs):
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
        vasp_input_set = vasp_input_set or MITMDSet(structure, start_temp=start_temp,
                                                    end_temp=end_temp, nsteps=nsteps,
                                                    **override_default_vasp_params)

        t = []
        if copy_vasp_outputs:
            t.append(CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"],
                                     contcar_to_poscar=True))

        t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, gamma_vasp_cmd=">>gamma_vasp_cmd<<",
                                  handler_group="md", wall_time=wall_time))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file,
                          additional_fields={"task_label": name}, defuse_unsuccessful=False))
        super(MDFW, self).__init__(t, parents=parents,
                                   name="{}-{}".format(structure.composition.reduced_formula, name),
                                   **kwargs)


class BoltztrapFW(Firework):
    def __init__(self, structure, name="boltztrap", db_file=None, parents=None, scissor=0.0,
                 soc=False, additional_fields=None, **kwargs):
        """
        Run Boltztrap (which includes writing bolztrap input files and parsing outputs). Assumes 
        you have a previous FW with the calc_locs passed into the current FW.

        Args:
            structure (Structure): - only used for setting name of FW
            name (str): name of this FW
            db_file (str): path to the db file
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            scissor (float): if scissor > 0, apply scissor on the band structure so that new
                band gap = scissor (in eV)
            soc (bool): whether the band structure is calculated with spin-orbit coupling
            additional_fields (dict): fields added to the document such as user-defined tags or name, ids, etc
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        additional_fields = additional_fields or {}
        t = [CopyVaspOutputs(calc_loc=True, contcar_to_poscar=True),
             RunBoltztrap(scissor=scissor, soc=soc),
             BoltztrapToDb(db_file=db_file, additional_fields=additional_fields),
             PassCalcLocs(name=name)]
        super(BoltztrapFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)


# TODO: @shyuep: Please do a code review before I look at this. -computron
class NEBRelaxationFW(Firework):
    """
    Relaxation Firework in NEB Workflow.

    Task 1) Read in a structure with "st_label" ("rlx", "ep0" or "ep1") and generates input sets.
    Task 2) Run VASP using Custodian
    Task 3) Update structure to spec
    Task 4) Pass CalcLocs named "{}_dir".format(st_label)
    """

    def __init__(self, spec, label, user_incar_settings=None, user_kpoints_settings=None,
                 additional_cust_args=None, **kwargs):
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
        if label == 'parent':
            vasp_input_set = MITRelaxSet(structure, user_incar_settings=user_incar_settings,
                                         user_kpoints_settings=user_kpoints_settings)
        else:  # label == "ep0" or "ep1"
            from pymatgen_diffusion.neb.io import MVLCINEBEndPointSet

            vasp_input_set = MVLCINEBEndPointSet(structure, user_incar_settings=user_incar_settings,
                                                 user_kpoints_settings=user_kpoints_settings)

        write_ep_task = WriteVaspFromIOSet(structure=structure, output_dir=".", vasp_input_set=vasp_input_set)

        # Task 2: Run VASP using Custodian
        cust_args = {"job_type": "normal", "gzip_output": False, "handler_group": "no_handler"}
        cust_args.update(additional_cust_args)
        run_vasp = RunVaspCustodian(vasp_cmd=">>vasp_cmd<<", gamma_vasp_cmd=">>gamma_vasp_cmd<<",
                                    **cust_args)

        # Task 3, 4: Transfer and PassCalLocs
        tasks = [write_ep_task, run_vasp, TransferNEBTask(label=label), PassCalcLocs(name=label)]

        super(NEBRelaxationFW, self).__init__(tasks, spec=spec, name=label, **kwargs)


# TODO: @shyuep: Please do a code review before I look at this. -computron
class NEBFW(Firework):
    """
    CI-NEB Firework in NEB Workflow.

    Task 1) Read in image structures from spec and generates input sets.
            The group of structures are labeled with neb_label (1, 2...)
    Task 2) Run NEB VASP using Custodian
    Task 3) Update structure to spec
    Task 4) Pass CalcLocs named "neb_{}".format(neb_label)
    """

    def __init__(self, spec, neb_label, from_images=True, user_incar_settings=None,
                 user_kpoints_settings=None, additional_cust_args=None, **kwargs):
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
            write_neb_task = WriteNEBFromImages(neb_label=neb_label,
                                                user_incar_settings=user_incar_settings,
                                                user_kpoints_settings=user_kpoints_settings)

        else:  # from endpoints
            write_neb_task = WriteNEBFromEndpoints(user_incar_settings=user_incar_settings,
                                                   user_kpoints_settings=user_kpoints_settings,
                                                   output_dir=".", sort_tol=sort_tol, d_img=d_img,
                                                   interpolation_type=interpolation_type)

        # Task 2: Run NEB using Custodian
        cust_args = {"job_type": "neb", "gzip_output": False, "handler_group": "no_handler"}
        cust_args.update(additional_cust_args)
        run_neb_task = RunVaspCustodian(vasp_cmd=">>vasp_cmd<<", gamma_vasp_cmd=">>gamma_vasp_cmd<<",
                                        **cust_args)

        # Task 3, 4: Transfer and PassCalcLocs
        tasks = [write_neb_task, run_neb_task, TransferNEBTask(label=label), PassCalcLocs(name=label)]

        super(NEBFW, self).__init__(tasks, spec=spec, name=label, **kwargs)
