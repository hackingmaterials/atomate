from typing import List, Optional, Union, Dict

import os

from atomate.common.firetasks.glue_tasks import (
    CopyFiles,
    CopyFilesFromCalcLoc,
    PassCalcLocs,
)

from atomate.vasp.firetasks.aimd import (
    CollectMDSegments,
    MDToDB
    )

from fireworks import Firework

__author__ = "Junsoo Park"
__email__ = "junsoo.park@nasa.gov"


class ARCMDFW(Firework):
    def __init__(
        self,
        structure,
        start_temp,
        end_temp,
        nsteps,
        name="AIMD",
        vasp_input_set=None,
        vasp_cmd=VASP_CMD,
        override_default_vasp_params=None,
        wall_time=19200,
        db_file=DB_FILE,
        parents=None,
        **kwargs,
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
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """
        override_default_vasp_params = override_default_vasp_params or {}
        vasp_input_set = vasp_input_set or ARCMDSet(
            structure,
            start_temp=start_temp,
            end_temp=end_temp,
            nsteps=nsteps,
            **override_default_vasp_params,
        )
        prev_runs = fw_spec.get("aimd", [])
        last_run = prev_runs[-1]
        last_config = last_run["final_configuration"]
        t = []
        t.append(WriteVaspFromIOSet(structure=last_config, vasp_input_set=vasp_input_set))
        t.append(
            RunVaspCustodian(
                vasp_cmd=vasp_cmd,
                gamma_vasp_cmd=">>vasp_cmd<<",#">>gamma_vasp_cmd<<",
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
        super().__init__(
            t,
            parents=parents,
            name=f"{structure.composition.reduced_formula}-{name}",
            **kwargs,
        )


class CollectMDSegmentsFW(Firework):
    """
    Compile serial AIMD runs segmented into multiple Fireworks

    Args:
        parents: Parent(s) of this Firework.
        name: Name of this FW.
        db_file: Path to a db file.
        cutoffs: A list of cutoffs to trial. If None, a set of trial cutoffs
            will be generated based on the structure (default).
        disp_cut: determines the mean displacement of perturbed
            structure to be included in harmonic (<) or anharmonic (>) fitting  
        bulk_modulus: in GPa, necessary for thermal expansion
        imaginary_tol: Tolerance used to decide if a phonon mode is
            imaginary, in THz.
        fit_method: Method used for fitting force constants. This can be
            any of the values allowed by the hiphive ``Optimizer`` class.
        mesh_density: The density of the q-point mesh used to calculate the
            phonon density of states. See the docstring for the ``mesh``
            argument in Phonopy.init_mesh() for more details.
        **kwargs: Other kwargs that are passed to Firework.__init__.
    """

    def __init__(
        self,
        name="Collect MD Runs",
        parents: Optional[Union[Firework, List[Firework]]] = None,
        db_file: str = None,
        **kwargs
    ):
        collect = CollectMDSegments()

        to_db = AIMDToDB(
            db_file=db_file, additional_fields={}
        )
        pass_locs = PassCalcLocs(name=name)

        tasks = [collect, to_db, pass_locs]
        super().__init__(tasks, parents=parents, name=name, **kwargs)


class FitForceConstantsFW(Firework):
    """
    Compile perturbed supercell calculations and fit force constants
    using hiPhive.

    Args:
        parents: Parent(s) of this Firework.
        name: Name of this FW.
        db_file: Path to a db file.
        cutoffs: A list of cutoffs to trial. If None, a set of trial cutoffs
            will be generated based on the structure (default).
        disp_cut: determines the mean displacement of perturbed
            structure to be included in harmonic (<) or anharmonic (>) fitting
        bulk_modulus: in GPa, necessary for thermal expansion
        imaginary_tol: Tolerance used to decide if a phonon mode is
            imaginary, in THz.
        fit_method: Method used for fitting force constants. This can be
            any of the values allowed by the hiphive ``Optimizer`` class.
        mesh_density: The density of the q-point mesh used to calculate the
            phonon density of states. See the docstring for the ``mesh``
            argument in Phonopy.init_mesh() for more details.
        **kwargs: Other kwargs that are passed to Firework.__init__.
    """

    def __init__(
            self,
            fit_method: str,
            disp_cut: float,
            bulk_modulus: float,
            imaginary_tol: float,
            mesh_density: float,
            temperature_qha: float,
            cutoffs: Optional[List[List[float]]] = None,
            name="Fit Force Constants",
            parents: Optional[Union[Firework, List[Firework]]] = None,
            db_file: str = None,
            **kwargs
    ):
        collect_structures = CollectPerturbedStructures()

        fit_force_constants = RunHiPhive(
            cutoffs=cutoffs,
            fit_method=fit_method,
            disp_cut=disp_cut,
            bulk_modulus=bulk_modulus,
            temperature_qha=temperature_qha,
            imaginary_tol=imaginary_tol,
        )
        to_db = ForceConstantsToDb(
            db_file=db_file, renormalized=False, mesh_density=mesh_density, additional_fields={}
        )
        pass_locs = PassCalcLocs(name=name)

        tasks = [collect_structures, fit_force_constants, to_db, pass_locs]
        super().__init__(tasks, parents=parents, name=name, **kwargs)
