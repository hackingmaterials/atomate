from fireworks import Firework
from pymatgen.io.vasp.sets import MPStaticSet

from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.vasp.config import DB_FILE, VASP_CMD
from atomate.vasp.firetasks import CopyVaspOutputs, RunVaspCustodian, VaspToDb
from atomate.vasp.firetasks.absorption_tasks import WriteVaspAbsorptionFromPrev
from atomate.vasp.firetasks.write_inputs import (
    WriteVaspFromIOSet,
    WriteVaspStaticFromPrev,
)


class AbsorptionFW(Firework):
    def __init__(
        self,
        structure,
        name="frequency dependent dielectrics",
        mode="STATIC",
        nbands=None,
        nbands_factor=2,
        reciprocal_density=200,
        nkred=None,
        nedos=2001,
        vasp_cmd=VASP_CMD,
        prev_calc_dir=None,
        db_file=DB_FILE,
        vasptodb_kwargs=None,
        parents=None,
        vasp_input_set_params=None,
        **kwargs,
    ):
        """
        FW that calculates frequency dependent dielectric function within
        indenpendent-particle-approxiamtion. A previous ground state calculation
        with the output WAVECAR is required by specifying mode = 'static'; in the case of no
        parent, a PBE functional ground state calculation will be performed and
        the WAVECAR will be saved. Then, perform another calculation with  'ALGO = EXACT, LOPTICS = True'
        with variable NBANDS by specifying MODE = "IPA". This calculation will save the
        WAVECAR and WAVEDER in case one wants to run RPA level absorption
        spectra. For RPA-DFT absorption spectrum, run another mode = 'RPA' calculation
        with the WAVECAR, WAVEDER saved from previous IPA calc.
        Args:
            structure (Structure): Input structure. For an interpolation, this
                is a dummy structure. See interpolate arg description.
            name (str): Name for the FireWork.
            mode: 'STATIC', 'IPA', or 'RPA'.
            nbands: number of bands to use, leave to None, and use nbands_factor instead
            nbands_factor: the multiplication of the number of bands
            reciprocal_density: k-point density
            nkred: reduced number of k-points, for RPA calculation use only, reduces the computing time
            nedos: energy mesh for DOS
            vasp_cmd (str): Command to run vasp.
            prev_calc_loc (bool or str): If true (default), copies outputs from previous calc. If
                a str value, retrieves a previous calculation output by name.
            vasp_input_set (str): string name for the VASP input set (e.g.,
                "MPAbsorptionSet").
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework. FW or list
                of FWS.
            vasp_input_set_params (dict): Dict of vasp_input_set_kwargs.
            prev_calc_dir (str): Path to a previous calculation to copy from
            vasptodb_kwargs (dict): kwargs to pass to VaspToDb
            **kwargs: Other kwargs that are passed to Firework.__init__.

        """
        t = []

        vasp_input_set_params = vasp_input_set_params or {}
        vasptodb_kwargs = vasptodb_kwargs or {}
        if "additional_fields" not in vasptodb_kwargs:
            vasptodb_kwargs["additional_fields"] = {}
        vasptodb_kwargs["additional_fields"]["task_label"] = name

        fw_name = "{}-{}-{}".format(
            structure.composition.reduced_formula if structure else "unknown",
            name,
            mode,
        )

        # define what wavecars to copy from the previous run
        if mode == "STATIC":
            wavecars = []
        elif mode == "IPA":
            wavecars = ["WAVECAR"]
        elif mode == "RPA":
            wavecars = ["WAVECAR", "WAVEDER"]
        else:
            raise Exception("Mode has to be from 'STATIC', 'IPA' or 'RPA'. ")

        #  "IPA" or "RPA" run
        if mode == "IPA" or mode == "RPA":
            if prev_calc_dir:
                # Copy the WAVECAR from previous calc directory
                t.append(
                    CopyVaspOutputs(
                        calc_dir=prev_calc_dir,
                        contcar_to_poscar=True,
                        additional_files=wavecars,
                    )
                )

                t.append(
                    WriteVaspAbsorptionFromPrev(
                        prev_calc_dir=".",
                        structure=structure,  # The structure will only be useful for the FW name
                        mode=mode,
                        copy_wavecar=True,
                        nbands=None,
                        nbands_factor=nbands_factor,
                        reciprocal_density=reciprocal_density,
                        nkred=nkred,
                        nedos=nedos,
                        **vasp_input_set_params,
                    )
                )

            elif parents:
                # Copy the WAVECAR from previous calc location
                t.append(
                    CopyVaspOutputs(
                        calc_loc=True, contcar_to_poscar=True, additional_files=wavecars
                    )
                )

                t.append(
                    WriteVaspAbsorptionFromPrev(
                        prev_calc_dir=".",
                        structure=structure,  # The structure will only be useful for the FW name
                        mode=mode,
                        copy_wavecar=True,
                        nbands=None,
                        nbands_factor=nbands_factor,
                        reciprocal_density=reciprocal_density,
                        nkred=nkred,
                        nedos=nedos,
                        **vasp_input_set_params,
                    )
                )

            else:
                raise ValueError(f"Must specify previous calculation for {mode}")

        # when mode = "static"
        elif mode == "STATIC":
            static_incar = {
                "LWAVE": True,
                "ENCUT": 500,
                "ISMEAR": 0,
                "SIGMA": 0.01,
                "LREAL": False,
                "GGA": "PE",
                "LELF": False,
                "LAECHG": False,
                "LASPH": False,
                "LVHAR": False,
                "LVTOT": False,
                "METAGGA": "None",
                "LMIXTAU": False,
            }

            if prev_calc_dir:
                # Copy only the CONTCAR from previous calc directory (often a relaxation run)
                t.append(
                    CopyVaspOutputs(
                        calc_dir=prev_calc_dir,
                        contcar_to_poscar=True,
                        additional_files=wavecars,
                    )
                )

                t.append(
                    WriteVaspStaticFromPrev(
                        reciprocal_density=reciprocal_density,
                        other_params={"user_incar_settings": static_incar},
                    )
                )

            elif parents:
                # Copy only the CONTCAR from previous calc
                t.append(
                    CopyVaspOutputs(
                        calc_loc=True, additional_files=wavecars, contcar_to_poscar=True
                    )
                )

                t.append(
                    WriteVaspStaticFromPrev(
                        reciprocal_density=reciprocal_density,
                        other_params={"user_incar_settings": static_incar},
                    )
                )

            elif structure:
                static_input_set = MPStaticSet(
                    structure=structure,
                    reciprocal_density=reciprocal_density,
                    user_incar_settings=static_incar,
                )

                t.append(
                    WriteVaspFromIOSet(
                        structure=structure, vasp_input_set=static_input_set
                    )
                )

            else:
                raise ValueError(
                    "Must specify structure or previous calculation for static calculation"
                )

        else:
            raise ValueError("Must specify a mode from 'STATIC', 'IPA', or 'RPA'")

        # use the 'default' custodian handler group
        handler_group = "default"

        # Run VASP
        t.append(
            RunVaspCustodian(
                vasp_cmd=vasp_cmd,
                auto_npar=">>auto_npar<<",
                handler_group=handler_group,
            )
        )
        t.append(PassCalcLocs(name=name))
        # Parse
        t.append(
            VaspToDb(
                db_file=db_file,
                additional_fields={
                    "task_label": structure.composition.reduced_formula
                    + " "
                    + name
                    + " "
                    + mode
                },
            )
        )

        super().__init__(t, parents=parents, name=fw_name, **kwargs)
