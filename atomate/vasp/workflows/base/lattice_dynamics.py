import math
import warnings
from copy import deepcopy
from typing import Dict, List, Optional, Union

import numpy as np

from atomate.utils.utils import get_logger
from atomate.vasp.config import DB_FILE, SHENGBTE_CMD, VASP_CMD
from atomate.vasp.firetasks import pass_vasp_result
from atomate.vasp.firetasks.lattice_dynamics import DEFAULT_TEMPERATURE
from atomate.vasp.fireworks.core import TransmuterFW
from atomate.vasp.fireworks.lattice_dynamics import (
    FitForceConstantsFW,
    LatticeThermalConductivityFW,
)
from atomate.vasp.powerups import add_additional_fields_to_taskdocs
from fireworks import Workflow
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.sets import MPStaticSet, VaspInputSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.advanced_transformations import (
    CubicSupercellTransformation,
)

__author__ = "Alex Ganose, Rees Chang"
__email__ = "aganose@lbl.gov, rc564@cornell.edu"
__date__ = "February 2020"

logger = get_logger(__name__)

_static_user_incar_settings = {
    "ADDGRID": True,
    "LCHARG": False,
    "ENCUT": 700,
    "ISMEAR": 0,
    "EDIFF": 1e-8,
    "PREC": "Accurate",
    "LAECHG": False,
    "LREAL": False,
    "LASPH": True,
}
_DEFAULT_SETTINGS = {"VASP_CMD": VASP_CMD, "DB_FILE": DB_FILE}
_WF_VERSION = 0.1


def get_lattice_dynamics_wf(
    structure: Structure,
    common_settings: Dict = None,
    vasp_input_set: Optional[VaspInputSet] = None,
    copy_vasp_outputs: bool = False,
    supercell_matrix_kwargs: Optional[dict] = None,
    num_supercell_kwargs: Optional[dict] = None,
    perturbed_structure_kwargs: Optional[dict] = None,
    calculate_lattice_thermal_conductivity: bool = False,
    thermal_conductivity_temperature: Union[float, Dict] = DEFAULT_TEMPERATURE,
    shengbte_cmd: str = SHENGBTE_CMD,
    shengbte_fworker: Optional[str] = None,
):
    """
    This workflow will calculate interatomic force constants and vibrational
    properties using the hiPhive package.

    A summary of the workflow is as follows:

    1. Calculate a supercell transformation matrix that brings the
       structure as close as cubic as possible, with all lattice lengths
       greater than 5 nearest neighbor distances.
    2. Perturb the atomic sites for each supercell using a Monte Carlo
       rattle procedure. The atoms are perturbed roughly according to a
       normal deviation. A number of standard deviation perturbation distances
       are included. Multiple supercells may be generated for each perturbation
       distance.
    3. Run static VASP calculations on each perturbed supercell to calculate
       atomic forces.
    4. Aggregate the forces and conduct the fit atomic force constants using
       the minimization schemes in hiPhive.
    5. Output the interatomic force constants, and phonon band structure and
       density of states to the database.
    6. Optional: Solve the lattice thermal conductivity using ShengBTE and
       output to the database.

    Args:
        structure: Initial structure.
        common_settings: Common settings dict. Supports "VASP_CMD", "DB_FILE",
            and "user_incar_settings" keys.
        vasp_input_set: Vasp input set for perturbed structure calculations.
        copy_vasp_outputs: Whether or not to copy previous vasp outputs.
        supercell_matrix_kwargs: Options that control the size of the supercell.
            Will be passed directly to CubicSupercellTransformation in
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

    if "n_configs_per_std" not in perturbed_structure_kwargs:
        n_supercells = get_num_supercells(supercell, **num_supercell_kwargs)
        perturbed_structure_kwargs["n_configs_per_std"] = n_supercells

    wf = get_perturbed_structure_wf(
        structure,
        supercell_matrix=supercell_matrix,
        vasp_input_set=vasp_input_set,
        common_settings=common_settings,
        copy_vasp_outputs=copy_vasp_outputs,
        pass_forces=True,
        **perturbed_structure_kwargs,
    )

    allow_fizzled = {"_allow_fizzled_parents": True}
    fw_fit_force_constant = FitForceConstantsFW(
        db_file=db_file, spec=allow_fizzled
    )
    wf.append_wf(Workflow.from_Firework(fw_fit_force_constant), wf.leaf_fw_ids)

    if calculate_lattice_thermal_conductivity:
        fw_lattice_conductivity = LatticeThermalConductivityFW(
            db_file=db_file,
            shengbte_cmd=shengbte_cmd,
            temperature=thermal_conductivity_temperature,
        )
        if shengbte_fworker:
            fw_lattice_conductivity.spec["_fworker"] = shengbte_fworker

        wf.append_wf(
            Workflow.from_Firework(fw_lattice_conductivity), [wf.fws[-1].fw_id]
        )

    formula = structure.composition.reduced_formula
    wf.name = "{} - lattice dynamics".format(formula)

    # Add workflow metadata to all relevant *ToDb tasks.
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
    copy_vasp_outputs: bool = False,
    rattle_stds: Optional[List[float]] = None,
    n_configs_per_std: int = 1,
    min_nn_scale: float = 0.85,
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
        copy_vasp_outputs: Whether or not to copy previous vasp outputs.
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
    override_vasp_params = {
        "user_incar_settings": common_settings["user_incar_settings"]
    }

    if supercell_matrix is None:
        supercell_matrix = np.eye(3)

    supercell_matrix = np.asarray(supercell_matrix).tolist()

    if rattle_stds is None:
        rattle_stds = np.linspace(0.01, 0.1, 5)

    if vasp_input_set is None:
        vasp_input_set = MPStaticSet(structure)

    # find the smallest nearest neighbor distance taking into account PBC
    min_distance = np.min(
        [d.nn_distance for d in structure.get_all_neighbors(10)[0]]
    )
    min_distance *= min_nn_scale
    all_rattle_stds = np.repeat(rattle_stds, n_configs_per_std)

    logger.debug("Using supercell_matrix of: {}".format(supercell_matrix))
    logger.debug(
        "Using {} supercells per displacement".format(n_configs_per_std)
    )
    logger.debug("Using {} rattle stds".format(len(rattle_stds)))
    logger.debug("Rattle stds: {}".format(rattle_stds))

    fws = []
    for i, rattle_std in enumerate(all_rattle_stds):
        fw_name = "{}: i={}; rattle_std={:.3f}".format(name, i, rattle_std)
        transformations = [
            "SupercellTransformation",
            "MonteCarloRattleTransformation",
        ]
        transformation_params = [
            {"scaling_matrix": supercell_matrix},
            {"rattle_std": rattle_std, "min_distance": min_distance},
        ]

        fw = TransmuterFW(
            name=fw_name,
            structure=structure,
            transformations=transformations,
            transformation_params=transformation_params,
            vasp_input_set=vasp_input_set,
            override_default_vasp_params=override_vasp_params,
            copy_vasp_outputs=copy_vasp_outputs,
            vasp_cmd=vasp_cmd,
            db_file=db_file,
        )

        if pass_forces:
            pass_dict = {
                "parent_structure": structure.to_json(),
                "supercell_matrix": supercell_matrix,
                "forces": ">>output.ionic_steps.-1.forces",
                "structure": ">>output.ionic_steps.-1.structure",
            }
            pass_task = pass_vasp_result(
                pass_dict=pass_dict,
                mod_spec_cmd="_push",
                mod_spec_key="perturbed_tasks",
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
    _, equiv_counts = np.unique(equiv_idxs, return_counts=True)
    min_count = np.min(equiv_counts)

    n_cells = math.ceil(min_num_equivalent_sites / min_count)
    return min([n_cells, max_num_supercells])


def _get_common_settings(common_settings: Optional[Dict]):
    common_settings = common_settings or {}
    for k, v in _DEFAULT_SETTINGS.items():
        if k not in common_settings:
            common_settings[k] = v

    user_incar_settings = deepcopy(_static_user_incar_settings)
    user_incar_settings.update(common_settings.get("user_incar_settings", {}))
    common_settings["user_incar_settings"] = user_incar_settings

    return common_settings
