from typing import List, Optional, Union

from atomate.common.firetasks.glue_tasks import (
    CopyFiles,
    CopyFilesFromCalcLoc,
    PassCalcLocs,
)
from atomate.vasp.analysis.lattice_dynamics import (
    FIT_METHOD,
    IMAGINARY_TOL,
    T_QHA,
    T_RENORM,
    T_KLAT
)
from atomate.vasp.config import SHENGBTE_CMD
from atomate.vasp.firetasks.lattice_dynamics import (
    CollectPerturbedStructures,
    ForceConstantsToDb,
    RunHiPhive,
    RunHiPhiveRenorm,
    RunShengBTE,
    ShengBTEToDb,
    MESH_DENSITY)
from fireworks import Firework

__author__ = "Alex Ganose, Junsoo Park"
__email__ = "aganose@lbl.gov, jsyony37@lbl.gov"


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
        separate_fit: If True, harmonic and anharmonic force constants are fit
            separately and sequentially, harmonic first then anharmonic. If
            False, then they are all fit in one go. Default is False.
        bulk_modulus: Must be supplied (in GPa) to copute thermal expansion
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
        name="Fit Force Constants",
        parents: Optional[Union[Firework, List[Firework]]] = None,
        db_file: str = None,
        cutoffs: Optional[List[List[float]]] = None,
        separate_fit: bool = False,
        bulk_modulus: float = None,
        imaginary_tol: float = IMAGINARY_TOL,
        fit_method: str = FIT_METHOD,
        mesh_density: float = MESH_DENSITY,
        **kwargs
    ):
        collect_structures = CollectPerturbedStructures()
        
        fit_force_constants = RunHiPhive(
            cutoffs=cutoffs,
            separate_fit=separate_fit,
            bulk_modulus=bulk_modulus,
            imaginary_tol=imaginary_tol,
            fit_method=fit_method
        )
        to_db = ForceConstantsToDb(
            db_file=db_file, renormalized=False, mesh_density=mesh_density, additional_fields={}
        )
        pass_locs = PassCalcLocs(name=name)

        tasks = [collect_structures, fit_force_constants, to_db, pass_locs]
        super().__init__(tasks, parents=parents, name=name, **kwargs)


class LatticeThermalConductivityFW(Firework):
    """
    Calculate the lattice thermal conductivity using ShengBTE.

    Args:
        name: Name of this FW.
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

    def __init__(
        self,
        name="Lattice Thermal Conductivity",
        prev_calc_dir: Optional[str] = None,
        db_file: str = None,
        shengbte_cmd: str = SHENGBTE_CMD,
        temperature: Union[float, dict] = T_KLAT,
        shengbte_control_kwargs: Optional[dict] = None,
        **kwargs
    ):
        # files needed to run ShengBTE
        files = [
            "structure_data.json",
            "FORCE_CONSTANTS_2ND",
            "FORCE_CONSTANTS_3RD",
        ]

        if prev_calc_dir:
            copy_files = CopyFiles(from_dir=prev_calc_dir, files_to_copy=files)
        else:
            copy_files = CopyFilesFromCalcLoc(calc_loc=True, filenames=files)

        run_shengbte = RunShengBTE(
            shengbte_cmd=shengbte_cmd,
            temperature=temperature,
            control_kwargs=shengbte_control_kwargs,
        )
        shengbte_to_db = ShengBTEToDb(db_file=db_file, additional_fields={})

        tasks = [copy_files, run_shengbte, shengbte_to_db]
        super().__init__(tasks, name=name, **kwargs)



class RenormalizationFW(Firework):        
    """
    Performs temperature-dependent phonon renormalization to obtain phonons
    at finite temperatures. Can be used to stabilize dynamically unstable modes
                                                                                                                                                                                                
    Args:
        name: Name of this FW.
        prev_calc_dir: Path to a directory containing the force constant
            information. Will override ``parents`` when collecting the force
            constants to run ShengBTE.
        db_file: Path to a db file.
        temperature: The temperature to perform phonon renormalization at
            Can be given as a single float, or a dictionary with the keys 
            "min", "max", "step".
        shengbte_control_kwargs: Options to be included in the ShengBTE
            control file.
        **kwargs: Other kwargs that are passed to Firework.__init__.
    """
    
    def __init__(
        self,
        name="Renormalization",
        prev_calc_dir: Optional[str] = None,
        db_file: str = None,
        temperature: Union[float, dict] = T_RENORM,
        **kwargs
    ):    
    
        # files needed to run renormalization
        files = ["cluster_space.cs","parameters.txt","force_constants.fcs"]
        
        if prev_calc_dir:
            copy_files = CopyFiles(from_dir=prev_calc_dir, files_to_copy=files)
        else:
            copy_files = CopyFilesFromCalcLoc(calc_loc=True, filenames=files)

        cs = ClusterSpace.read('cluster_space.cs')
        fcs = ForceConstants.read('force_constants.fcs')
        param = np.loadtxt('parameters.txt')
        
        renorm_force_constants = RunHiPhiveRenorm(
            cs=cs,
            param=param,
            fcs=fcs,
            T_renorm=T_renorm
        )        

        to_db = ForceConstantsToDb(
            db_file=db_file, renormalized=True, mesh_density=mesh_density, additional_fields={}
	)
        pass_locs = PassCalcLocs(name=name)

        tasks = [renorm_force_constants, to_db, pass_locs]
        super().__init__(tasks, name=name, **kwargs)
