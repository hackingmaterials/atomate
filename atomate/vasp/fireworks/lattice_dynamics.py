from typing import Optional, List, Union

from atomate.common.firetasks.glue_tasks import PassCalcLocs, \
    CopyFilesFromCalcLoc, CopyFiles
from atomate.vasp.analysis.lattice_dynamics import IMAGINARY_TOL, \
    MAX_N_IMAGINARY, MAX_IMAGINARY_FREQ, FIT_METHOD
from atomate.vasp.analysis.phonopy import MESH_DENSITY
from atomate.vasp.config import SHENGBTE_CMD
from atomate.vasp.firetasks.lattice_dynamics import CollectPerturbedStructures, \
    RunHiPhive, ForceConstantsToDb, DEFAULT_TEMPERATURE, ShengBTEToDb, \
    RunShengBTE
from fireworks import Firework

__author__ = "Alex Ganose"
__email__ = "aganose@lbl.gov"


class FitForceConstantsFW(Firework):

    def __init__(
        self,
        name="Fit Force Constants",
        parents: Optional[Union[Firework, List[Firework]]] = None,
        db_file: str = None,
        cutoffs: Optional[List[List[float]]] = None,
        imaginary_tol: float = IMAGINARY_TOL,
        max_n_imaginary: int = MAX_N_IMAGINARY,
        max_imaginary_freq: float = MAX_IMAGINARY_FREQ,
        fit_method: str = FIT_METHOD,
        mesh_density: float = MESH_DENSITY,
        **kwargs
    ):
        """
        Compile perturbed supercell calculations and fit force constants
        using hiPhive.

        Args:
            parents: Parent(s) of this Firework.
            name: Name of this FW.
            db_file: Path to a db file.
            cutoffs: A list of cutoffs to trial. If None, a set of trial cutoffs
                will be generated based on the structure (default).
            imaginary_tol: Tolerance used to decide if a phonon mode is
                imaginary, in THz.
            max_n_imaginary: Maximum number of imaginary modes allowed in the
                final fitted force constant solution. If this criteria is not
                reached by any cutoff combination the Firework will fizzle.
            max_imaginary_freq: Maximum allowed imaginary frequency in the
                final fitted force constant solution. If this criteria is not
                reached by any cutoff combination this FireTask will fizzle.
            fit_method: Method used for fitting force constants. This can be
                any of the values allowed by the hiphive ``Optimizer`` class.
            mesh_density: The density of the q-point mesh used to calculate the
                phonon density of states. See the docstring for the ``mesh``
                argument in Phonopy.init_mesh() for more details.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """
        collect_structures = CollectPerturbedStructures(db_file)
        fit_constants = RunHiPhive(
            cutoffs=cutoffs,
            imaginary_tol=imaginary_tol,
            max_n_imaginary=max_n_imaginary,
            max_imaginary_freq=max_imaginary_freq,
            fit_method=fit_method
        )
        to_db = ForceConstantsToDb(db_file=db_file, mesh_density=mesh_density)
        pass_locs = PassCalcLocs(name=name)

        tasks = [collect_structures, fit_constants, to_db, pass_locs]
        super().__init__(tasks, parents=parents, name=name, **kwargs)


class LatticeThermalConductivityFW(Firework):

    def __init__(
            self,
            name="Lattice Thermal Conductivity",
            parents: Optional[Union[Firework, List[Firework]]] = None,
            prev_calc_dir: Optional[str] = None,
            db_file: str = None,
            shengbte_cmd: str = SHENGBTE_CMD,
            temperature: Union[float, dict] = DEFAULT_TEMPERATURE,
            shengbte_control_kwargs: Optional[dict] = None,
            **kwargs
    ):
        """
        Calculate the lattice thermal conductivity using ShengBTE.

        Args:
            name: Name of this FW.
            parents: Parent(s) of this Firework.
            prev_calc_dir: Path to a directory containing the force constant
                information. Will override ``parents`` when collecting the force
                constants to run ShengBTE.
            db_file: Path to a db file.
            shengbte_cmd: The name of the shengbte executable to run. Supports
                env_chk.
            temperature: The temperature to calculate the lattice thermal
                conductivity for. Can be given as a single float, or a
                dictionary with the keys "min", "max", "step".
            shengbte_control_kwargs: Options to be included in the ShengBTE
                control file.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """
        # files needed to run ShengBTE
        files = [
            "structure_data.json", "FORCE_CONSTANTS_2ND", "FORCE_CONSTANTS_3RD"
        ]

        if prev_calc_dir:
            copy_files = CopyFiles(from_dir=prev_calc_dir, files_to_copy=files)
        elif parents:
            copy_files = CopyFilesFromCalcLoc(calc_loc=True, filenames=files)
        else:
            raise ValueError("Must specify parents or prev_calc_dir.")

        run_shengbte = RunShengBTE(
            shengbte_cmd=shengbte_cmd,
            temperature=temperature,
            control_kwargs=shengbte_control_kwargs
        )
        shengbte_to_db = ShengBTEToDb(db_file=db_file)

        tasks = [copy_files, run_shengbte, shengbte_to_db]
        super().__init__(tasks, parents=parents, name=name, **kwargs)
