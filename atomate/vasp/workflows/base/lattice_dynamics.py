import math

import numpy as np

from copy import deepcopy

from monty.dev import requires
from typing import List, Dict, Optional, Union

from uuid import uuid4

from atomate.vasp.firetasks.lattice_dynamics import FitForceConstants, \
    ForceConstantsToDB, RunShengBTE, ShengBTEToDB
from pymatgen.io.vasp.sets import MPStaticSet
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.advanced_transformations import (
    CubicSupercellTransformation,
)

from fireworks import Workflow, Firework

from atomate.utils.utils import get_logger
from atomate.vasp.config import VASP_CMD, DB_FILE, SHENGBTE_CMD
from atomate.vasp.fireworks.core import StaticFW
from atomate.vasp.powerups import add_additional_fields_to_taskdocs
from atomate.vasp.analysis.lattice_dynamics import generate_perturbed_structures

try:
    import hiphive
except ImportError:
    hiphive = False

__author__ = "Rees Chang, Alex Ganose"
__email__ = "rc564@cornell.edu, aganose@lbl.gov"
__date__ = "July 2019"
__lattice_dynamics_wf_version__ = 1.0

logger = get_logger(__name__)

_static_uis = {
    "ADDGRID": True,  # Fast Fourier Transform grid
    "LCHARG": False,
    "ENCUT": 700,
    "EDIFF": 1e-8,  # may need to tune this
    "PREC": "Accurate",
    "LAECHG": False,
    "LREAL": False,
    "LASPH": True,
}

# define shared constants
MIN_NN_SCALE = 0.85
MAX_IMAGINARY_FREQ = 10  # in THz
IMAGINARY_TOL = 0.025  # in THz
MAX_N_IMAGINARY = 3
DEFAULT_TEMPERATURE = 300
FIT_METHOD = "least-squares"
MESH_DENSITY = 100.  # should always be a float
_DEFAULT_SETTINGS = {"VASP_CMD": VASP_CMD, "DB_FILE": DB_FILE}


@requires(hiphive, "hiphive is required for lattice dynamics workflow")
def get_lattice_dynamics_wf(
    structure: Structure,
    common_settings: Dict = None,
    symmetrize: bool = False,
    symprec: float = 0.01,
    min_atoms: Optional[int] = None,
    max_atoms: Optional[int] = None,
    num_nn_dists: float = 5.,
    force_diagonal_transformation: bool = True,
    rattle_stds: Optional[List[float]] = None,
    min_nn_scale: float = MIN_NN_SCALE,
    min_num_equivalent_sites: int = 100,
    max_num_supercells: int = 30,
    dynamic_static_calcs: bool = False,
    calculate_lattice_thermal_conductivity: bool = False,
    thermal_conductivity_temperature: Union[float, Dict] = 300,
    shengbte_cmd: str = SHENGBTE_CMD,
    shengbte_fworker: Optional[str] = None,
):
    """
    This workflow will calculate interatomic force constants and vibrational
    properties using the hiPhive package.

    A summary of the workflow is as follows:

    1. Transform the input structure into a supercell
    2. Transform the supercell into a list of supercells with all atoms
       randomly perturbed from their original sites
    3. Run static VASP calculations on each perturbed supercell to
       calculate atomic forces.
    4. Aggregate the forces and conduct the fit the force constants using
       the minimization schemes in hiPhive.
    5. Output the interatomic force constants, and phonon band structure and
       density of states to the database.
    6. Optional: Solve the lattice thermal conductivity using ShengBTE.

    Args:
        structure: Initial structure.
        common_settings: Common settings dict.
        symmetrize: Whether to symmetrize the structure before running the
            workflow.
        symprec: Symmetry precision for determining the symmetry of the
            structure.
        min_atoms: Minimum number of atoms to constrain the supercell.
        max_atoms: Maximum number of atoms to constrain the supercell.
        num_nn_dists: Number of nearest neighbor distances that the shortest
            direction of the supercell must be as long as.
        force_diagonal_transformation: If True, the supercell
            transformation will be constrained to have a diagonal
            transformation matrix. If False, the supercell transformation
            will not be constrained to be diagonal (resulting in a more
            cubic supercell). A diagonal supercell is required to calculate
            lattice thermal conductivity.
        rattle_stds: List of standard deviations (in normal distribution) to
            to use when displacing the sites.
        min_nn_scale: Controls the minimum interatomic distance used for
            computing the probability for each rattle move.
        min_num_equivalent_sites: The minimum number of equivalent sites for
            each species. The aim is to have each equivalent site represented at
            least ``min_num_equivalent_sites`` times. For example, if a symmetry
            inequivalent site appears 26 times in the supercell structure,
            and ``min_num_equivalent_sites`` is 100, then at least 4 supercells
            are required for that site to appear over 100 times.
        max_num_supercells: The maximum number of supercells per displacement.
            For very unsymmetric structures, this imposes a hard limit on
            the number of sites generated due to ``min_num_equivalent_sites``.
        dynamic_static_calcs: If True and hiPhive fails to return all real
             harmonic phonon frequencies on the first try, then the workflow
             will dynamically create more static perturbed supercell
             calculations with larger displacements to improve the fitting.
        calculate_lattice_thermal_conductivity: If True and force constant
            fitting does not return imaginary modes, then use ShengBTE to
            calculate the lattice thermal conductivity.
        thermal_conductivity_temperature: The temperature to calculate the
            lattice thermal conductivity for. Can be given as a single float, or
            a dictionary with the keys "min", "max", "step".
        shengbte_cmd: Command to run ShengBTE. Supports env_chk.
        shengbte_fworker: If None, the ShengBTE firework's fworker will be set
            to all the previous fireworks' fworker. If str, the ShengBTE
            firework's fworker will be set to shengbte_fworker.
   """
    common_settings = _get_common_settings(common_settings)
    uis = deepcopy(_static_uis)
    uis.update(common_settings.get("user_incar_settings", {}))

    if not force_diagonal_transformation and calculate_lattice_thermal_conductivity:
        raise ValueError(
            "Diagonal transformation required to calculate lattice thermal "
            "conductivity"
        )

    sga = SpacegroupAnalyzer(structure, symprec=symprec)
    if symmetrize:
        structure = sga.get_primitive_standard_structure()

    st = CubicSupercellTransformation(
        min_atoms,
        max_atoms,
        num_nn_dists,
        force_diagonal_transformation=force_diagonal_transformation,
    )
    supercell = st.apply_transformation(structure)

    n_cells = _get_n_cells(sga, min_num_equivalent_sites, max_num_supercells)
    perturbed_supercells, disp_dists, fws = get_perturbed_structure_fws(
        supercell,
        common_settings=common_settings,
        user_incar_settings=uis,
        rattle_stds=rattle_stds,
        min_nn_scale=min_nn_scale,
        n_configs_per_std=n_cells,
    )

    logger.debug("Generating {} supercells per displacement".format(n_cells))
    logger.debug("Using {} displacements".format(len(disp_dists)))
    logger.debug("Displacements: {}".format(disp_dists))

    wf_meta = {
        "wf_uuid": str(uuid4()),
        "wf_name": "LatticeDynamicsWF",
        "wf_version": __lattice_dynamics_wf_version__,
    }

    # Add tasks to fit force constants
    fit_forces = FitForceConstants(
        db_file=common_settings["DB_FILE"],
        wf_uuid=wf_meta["wf_uuid"],
        parent_structure=structure,
        supercell_structure=supercell,
        supercell_matrix=st.transformation_matrix.tolist(),
    )
    forces_to_db = ForceConstantsToDB(
        db_file=common_settings["DB_FILE"], wf_uuid=wf_meta["wf_uuid"],
    )
    fws.append(
        Firework(
            [fit_forces, forces_to_db],
            name="Force Constant Fitting",
            parents=fws[-len(perturbed_supercells):],
        )
    )

    # Add ShengBTE tasks
    run_shengbte = RunShengBTE(
        db_file=common_settings["DB_FILE"],
        wf_uuid=wf_meta["wf_uuid"],
        structure=structure,
        supercell_matrix=st.transformation_matrix.tolist(),
        shengbte_cmd=shengbte_cmd,
        temperature=thermal_conductivity_temperature,
    )
    shengbte_to_db = ShengBTEToDB(
        db_file=common_settings["DB_FILE"],
        wf_uuid=wf_meta["wf_uuid"],
    )
    fws.append(
        Firework(
            [run_shengbte, shengbte_to_db],
            name="Lattice Thermal Conductivity",
            parents=fws[-1]
        )
    )

    formula = structure.composition.reduced_formula
    wf_name = "{} - lattice dynamics".format(formula)

    wf = Workflow(fws, name=wf_name)
    wf = add_additional_fields_to_taskdocs(
        wf, {"wf_meta": wf_meta}, task_name_constraint="VaspToDb"
    )

    return wf


def get_perturbed_structure_fws(
    supercell: Structure,
    common_settings: Optional[Dict] = None,
    user_incar_settings: Optional[Dict] = None,
    rattle_stds: Optional[List[float]] = None,
    n_configs_per_std: int = 1,
    min_nn_scale: float = MIN_NN_SCALE,
) -> List[Firework]:
    """
    Get static calculations to calculate the forces on each perturbed supercell.

    Args:
        supercell: Parent supercell structure.
        common_settings: Common settings dict
        user_incar_settings: User incar settings override.
        rattle_stds: List of standard deviations (in normal distribution) to
            to use when displacing the sites.
        n_configs_per_std: Number of structures to generate per rattle standard
            deviation.
        min_nn_scale: Controls the minimum interatomic distance used for
            computing the probability for each rattle move.

    Returns:
        A list of static fireworks.
    """
    user_incar_settings = user_incar_settings if user_incar_settings else {}
    common_settings = _get_common_settings(common_settings)

    # Generate list of perturbed supercells
    perturbed_supercells, disp_dists = generate_perturbed_structures(
        supercell,
        rattle_stds=rattle_stds,
        n_configs_per_std=n_configs_per_std,
        min_nn_scale=min_nn_scale
    )

    fws = []
    data = enumerate(zip(perturbed_supercells, disp_dists))
    for i, (supercell, dist) in data:
        name = "perturbed supercell, idx: {}, disp_val: {:.3f} static".format(
            i, dist
        )

        vis = MPStaticSet(supercell, user_incar_settings=user_incar_settings)
        static_fw = StaticFW(
            supercell,
            vasp_input_set=vis,
            vasp_cmd=common_settings["VASP_CMD"],
            db_file=common_settings["DB_FILE"],
            name=name,
        )
        static_fw.spec["displacement_value"] = dist
        fws.append(static_fw)

    return fws


def _get_common_settings(common_settings: Optional[Dict]):
    common_settings = common_settings or {}
    for k, v in _DEFAULT_SETTINGS.items():
        if k not in common_settings:
            common_settings[k] = v
    return common_settings


def _get_n_cells(
    spacegroup_analyzer: SpacegroupAnalyzer,
    min_num_equivalent_sites: int,
    max_num_structure: int,
) -> int:
    """
    Get the number of supercells to run per rattle standard deviation.

    The aim is to have each equivalent site represented at least
    ``min_num_equivalent_sites`` times. For example, if a symmetry inequivalent
    site appears 26 times in the supercell structure, `n_configs` would be 4.
    I.e., 4 * 26 > 100.

    Args:
        spacegroup_analyzer: A spacegroup analyzer object.
        min_num_equivalent_sites: The minimum number of equivalent sites.
        max_num_structure: The maximum number of cells.

    Returns:
        The number of supercells needed for all symmetry inequivalent sites to
        be represented at least `min_num_equivalent_sites` times.
    """
    # get the equivalent sites in the structure and calculate the smallest
    # site degeneracy
    equiv_idxs = spacegroup_analyzer.get_symmetry_dataset()["equivalent_atoms"]
    equiv_counts = np.unique(equiv_idxs, return_counts=True)
    min_count = np.min(equiv_counts)

    n_cells = math.ceil(min_num_equivalent_sites / min_count)
    return min(n_cells, max_num_structure)
