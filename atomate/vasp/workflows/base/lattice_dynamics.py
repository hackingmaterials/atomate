import warnings
import math

import numpy as np

from copy import deepcopy
from typing import List, Dict, Optional, Union


from pymatgen.io.vasp.sets import MPStaticSet, VaspInputSet
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.advanced_transformations import (
    CubicSupercellTransformation,
)

from monty.dev import requires
from fireworks import Workflow

from atomate.utils.utils import get_logger
from atomate.vasp.config import VASP_CMD, DB_FILE, SHENGBTE_CMD
from atomate.vasp.fireworks.core import TransmuterFW
from atomate.vasp.powerups import add_additional_fields_to_taskdocs
from atomate.vasp.firetasks import pass_vasp_result
from atomate.vasp.fireworks.lattice_dynamics import FitForceConstantsFW, \
    LatticeThermalConductivityFW

try:
    import hiphive
except ImportError:
    hiphive = False

__author__ = "Rees Chang, Alex Ganose"
__email__ = "rc564@cornell.edu, aganose@lbl.gov"
__date__ = "July 2019"

logger = get_logger(__name__)

_static_user_incar_settings = {
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
_WF_VERSION = 1.0


@requires(hiphive, "hiphive is required for lattice dynamics workflow")
def get_lattice_dynamics_wf(
    structure: Structure,
    common_settings: Dict = None,
    vasp_input_set: Optional[VaspInputSet] = None,
    supercell_matrix_kwargs: Optional[dict] = None,
    num_supercell_kwargs: Optional[dict] = None,
    perturbed_structure_kwargs: Optional[dict] = None,
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
        common_settings: Common settings dict. Supports "VASP_CMD", "DB_FILE",
            and "user_incar_settings" keys.
        vasp_input_set: Vasp input set for perturbed structure calculations.
        supercell_matrix_kwargs: Options that control the size of the supercell
            that will be perturbed. Will be passed directly to
            CubicSupercellTransformation in
            pymatgen.transformations.advanced_transformations. Note, a diagonal
            supercell is required to calculate lattice thermal conductivity.
        num_supercell_kwargs: Options that control the number of supercells
            generated per perturbation standard deviation distance. See the
            docstring of ``get_num_supercells`` for more details.
        perturbed_structure_kwargs: Options that control the number of
            and magnitude of atomic displacements. Currently, the options are
            "rattle_std", "n_configs_per_std", and "min_nn_scale". See the
            docstring for ``get_perturbed_structure_wf`` for more details.
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
    supercell_matrix_kwargs = supercell_matrix_kwargs or {}
    num_supercell_kwargs = num_supercell_kwargs or {}
    perturbed_structure_kwargs = perturbed_structure_kwargs or {}
    common_settings = _get_common_settings(common_settings)
    db_file = common_settings["DB_FILE"]

    if calculate_lattice_thermal_conductivity:
        if supercell_matrix_kwargs.get("force_diagonal", True):
            warnings.warn(
                "Diagonal transformation required to calculate lattice thermal "
                "conductivity. Forcing diagonal matrix"
            )
        supercell_matrix_kwargs["force_diagonal"] = True

    st = CubicSupercellTransformation(**supercell_matrix_kwargs)
    supercell = st.apply_transformation(structure)
    supercell_matrix = st.transformation_matrix

    if "n_supercells" not in perturbed_structure_kwargs:
        n_supercells = get_num_supercells(supercell, **num_supercell_kwargs)
        perturbed_structure_kwargs["n_supercells"] = n_supercells

    wf = get_perturbed_structure_wf(
        structure,
        supercell_matrix=supercell_matrix,
        name="ld perturbed structure",
        vasp_input_set=vasp_input_set,
        common_settings=common_settings,
        pass_forces=True,
        **perturbed_structure_kwargs,
    )

    allow_fizzled = {"_allow_fizzled_parents": True}
    fw_fit_force_constant = FitForceConstantsFW(
        parents=wf.fws[wf.leaf_fw_ids], db_file=db_file, spec=allow_fizzled
    )
    wf.fws.append(fw_fit_force_constant)

    if calculate_lattice_thermal_conductivity:
        fw_lattice_conductivity = LatticeThermalConductivityFW(
            parents=wf.fws[-1], db_file=db_file
        )
        wf.fws.append(fw_lattice_conductivity)

    formula = structure.composition.reduced_formula
    wf.name = "{} - lattice dynamics".format(formula)

    # Add workflow meta data to all relevant *ToDb tasks.
    wf_meta = {"wf_name": "LatticeDynamicsWF", "wf_version": _WF_VERSION}
    for task_name in ["VaspToDb", "ShengBTEToDb", "ForceConstantsToDb"]:
        wf = add_additional_fields_to_taskdocs(
            wf, {"wf_meta": wf_meta}, task_name_constraint=task_name
        )

    return wf


def get_perturbed_structure_wf(
    structure: Structure,
    supercell_matrix: Optional[np.ndarray],
    name: str = "perturbed structure",
    vasp_input_set: Optional[VaspInputSet] = None,
    common_settings: Optional[Dict] = None,
    rattle_stds: Optional[List[float]] = None,
    n_configs_per_std: int = 1,
    min_nn_scale: float = MIN_NN_SCALE,
    pass_forces: bool = True,
) -> Workflow:
    """
    Get static calculations to calculate the forces on each perturbed supercell.

    Args:
        structure: Input structure.
        supercell_matrix: Supercell transformation matrix
        name: Transmuter firework name.
        vasp_input_set: Vasp input set for perturbed structure calculations.
        common_settings: Common settings dict. Supports "VASP_CMD", "DB_FILE",
            and "user_incar_settings" keys.
        rattle_stds: List of standard deviations (in normal distribution) to
            to use when displacing the sites.
        n_configs_per_std: Number of structures to generate per rattle standard
            deviation.
        min_nn_scale: Controls the minimum interatomic distance used for
            computing the probability for each rattle move.
        pass_forces: Whether to append the force and supercell structure
            information into the perturbed_tasks key of the child fireworks.

    Returns:
        The workflow.
    """
    common_settings = _get_common_settings(common_settings)
    vasp_cmd = common_settings["VASP_CMD"]
    db_file = common_settings["DB_FILE"]
    user_incar_settings = common_settings["user_incar_settings"]

    if supercell_matrix is None:
        supercell_matrix = np.eye(3)

    supercell_matrix = np.asarray(supercell_matrix).tolist()

    if rattle_stds is None:
        rattle_stds = np.linspace(0.01, 0.1, 5)

    if vasp_input_set is None:
        vasp_input_set = MPStaticSet(structure)
    else:
        # ensure we don't override the user_incar_settings in the input set
        user_incar_settings.update(vasp_input_set.user_incar_settings)
    vasp_input_set.user_incar_settings = user_incar_settings

    min_distance = np.min(structure.distance_matrix) * min_nn_scale
    all_rattle_stds = np.repeat(rattle_stds, n_configs_per_std)

    logger.debug("Using {} supercells / displacement".format(n_configs_per_std))
    logger.debug("Using {} rattle stds".format(len(rattle_stds)))
    logger.debug("Rattle stds: {}".format(rattle_stds))

    fws = []
    for i, rattle_std in all_rattle_stds:
        name = "{} : i = {}; rattle_std : {:.3f}".format(name, i, rattle_std)
        transformations = [
            "SupercellTransformation", "MonteCarloRattleTransformation"
        ]
        transformation_params = [
            {"scaling_matrix": supercell_matrix},
            {"rattle_std": rattle_std, "min_distance": min_distance}
        ]

        fw = TransmuterFW(
            name=name,
            structure=structure,
            transformations=transformations,
            transformation_params=transformation_params,
            vasp_input_set=vasp_input_set,
            copy_vasp_outputs=True,
            vasp_cmd=vasp_cmd,
            db_file=db_file,
        )

        if pass_forces:
            pass_dict = {
                'parent_structure': structure.to_json(),
                'supercell_matrix': supercell_matrix,
                'forces': '>>output.ionic_steps.-1.forces',
                'structure': '>>output.ionic_steps.-1.structure',
            }
            pass_task = pass_vasp_result(
                pass_dict=pass_dict,
                mod_spec_cmd="_push",
                mod_spec_key="perturbed_tasks"
            )
            fw.tasks.append(pass_task)

        fws.append(fw)

    wfname = "{}: {}".format(structure.composition.reduced_formula, name)
    return Workflow(fws, name=wfname)


def get_num_supercells(
    supercell_structure: Structure,
    symprec: float = 0.01,
    min_num_equivalent_sites: int = 100,
    max_num_supercells: int = 30,
) -> int:
    """
    Get the number of supercells to run per rattle standard deviation.

    The aim is to have each equivalent site represented at least
    ``min_num_equivalent_sites`` times. For example, if a symmetry inequivalent
    site appears 26 times in the supercell structure, `n_configs` would be 4.
    I.e., 4 * 26 > 100.

    Args:
        supercell_structure: The ideal supercell structure.
        symprec: The symmetry precision for determining if sites are equivalent.
        min_num_equivalent_sites: The minimum number of equivalent sites.
        max_num_supercells: The maximum number of cells.

    Returns:
        The number of supercells needed for all symmetry inequivalent sites to
        be represented at least `min_num_equivalent_sites` times.
    """
    # get the equivalent sites in the structure and calculate the smallest
    # site degeneracy
    sga = SpacegroupAnalyzer(supercell_structure, symprec=symprec)
    equiv_idxs = sga.get_symmetry_dataset()["equivalent_atoms"]
    equiv_counts = np.unique(equiv_idxs, return_counts=True)
    min_count = np.min(equiv_counts)

    n_cells = math.ceil(min_num_equivalent_sites / min_count)
    return min(n_cells, max_num_supercells)


def _get_common_settings(common_settings: Optional[Dict]):
    common_settings = common_settings or {}
    for k, v in _DEFAULT_SETTINGS.items():
        if k not in common_settings:
            common_settings[k] = v

    user_incar_settings = deepcopy(_static_user_incar_settings)
    user_incar_settings.update(common_settings.get("user_incar_settings", {}))
    common_settings["user_incar_settings"] = user_incar_settings

    return common_settings
