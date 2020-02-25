# coding: utf-8
from copy import deepcopy
from typing import List, Dict, Optional

from uuid import uuid4

from pymatgen.transformations.advanced_transformations import (
    CubicSupercellTransformation
)
from pymatgen.io.vasp.sets import MPStaticSet
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from fireworks import Workflow, Firework
from atomate.utils.utils import get_logger
from atomate.vasp.config import VASP_CMD, DB_FILE
from atomate.vasp.fireworks.core import StaticFW
from atomate.vasp.firetasks.parse_outputs import CSLDForceConstantsToDB
from atomate.vasp.powerups import add_additional_fields_to_taskdocs
from atomate.vasp.analysis.csld import generate_perturbed_supercells

__author__ = "Rees Chang, Alex Ganose"
__email__ = "rc564@cornell.edu, aganose@lbl.gov"
__date__ = "July 2019"
__csld_wf_version__ = 1.0

logger = get_logger(__name__)

_static_uis = {
    "ADDGRID": True,  # Fast Fourier Transform grid
    "LCHARG": False,
    "ENCUT": 700,
    "EDIFF": 1e-8,  # may need to tune this
    "PREC": 'Accurate',
    "LAECHG": False,
    "LREAL": False,
    "LASPH": True
}

# define shared constants
_MAX_DISPLACEMENT = 0.1
_MIN_DISPLACEMENT = 0.01
_NUM_DISPLACEMENTS = 10
_MIN_RANDOM_DISPLACEMENTS = None
_CELLS_PER_DISPLACEMENT = 1
_DEFAULT_SETTINGS = {"VASP_CMD": VASP_CMD, "DB_FILE": DB_FILE}


def get_lattice_dynamics_wf(
    structure: Structure,
    common_settings: Dict = None,
    symprec: Optional[float] = None,
    min_atoms: Optional[int] = None,
    max_atoms: Optional[int] = None,
    num_nn_dists: float = 5,
    force_diagonal_transformation: bool = True,
    max_displacement: float = _MAX_DISPLACEMENT,
    min_displacement: float = _MIN_DISPLACEMENT,
    num_displacements: int = _NUM_DISPLACEMENTS,
    n_cells_per_displacement_distance: int = _CELLS_PER_DISPLACEMENT,
    min_random_distance: Optional[float] = _MIN_RANDOM_DISPLACEMENTS,
    dynamic_static_calcs: bool = False,
    do_shengbte: bool = False,
    shengbte_t_range: Optional[List[float]] = None,
    shengbte_fworker: Optional[str] = None
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
        symprec: Symmetry precision for symmetrizing the structure. If set to
            ``None`` (the default), the structure will not be symmetrized.
        min_atoms: Minimum number of atoms to constrain the supercell.
        max_atoms: Maximum number of atoms to constrain the supercell.
        num_nn_dists: Number of nearest neighbor distances that the shortest
            direction of the supercell must be as long as.
        force_diagonal_transformation: If True, the supercell
            transformation will be constrained to have a diagonal
            transformation matrix. If False, the supercell transformation
            will not be constrained to be diagonal (resulting in a more
            cubic supercell).
        max_displacement: Maximum displacement distance for perturbing the
            sites in Angstroms.
        min_displacement: Minimum displacement distance for perturbing the
            sites in Angstroms.
        num_displacements: Number of unique displacement distances to
            try, uniformly distributed between ``min_displacement`` and
            ``max_displacement``.
        n_cells_per_displacement_distance: Number of perturbed
            structures to generate for each unique displacement distance.
        min_random_distance: If None (default), then for a given perturbed
            structure, all atoms will move the same distance from their original
            locations. If float, then for a given perturbed structure, the
            distances the atoms move will be uniformly distributed from a
            minimum distance of ``min_random_distance`` to one of the
            displacement distances uniformly sampled between
            ``min_displacement`` and ``max_displacement``.
        dynamic_static_calcs: If True and hiPhive fails to return all real
             harmonic phonon frequencies on the first try, then the workflow
             will dynamically create more static perturbed supercell
             calculations with larger displacements to improve the fitting.
        do_shengbte: If True and force constant fitting does not return
            imaginary modes, then the workflow will dynamically create a
            ShengBTE firework for calculating the lattice thermal conductivity.
        shengbte_t_range: A list of temperatures at which to calculate the
            lattice thermal conductivity using ShengBTE. If None (default),
            then 300 K will be used.
        shengbte_fworker: If None, the ShengBTE firework's fworker will be set
            to all the previous fireworks' fworker. If str, the ShengBTE
            firework's fworker will be set to shengbte_fworker.
   """
    common_settings = common_settings or _DEFAULT_SETTINGS
    uis = deepcopy(_static_uis)
    uis.update(common_settings.get("user_incar_settings", {}))

    if force_diagonal_transformation is True and num_nn_dists == 5:
        num_nn_dists = 6

    if symprec is not None:
        sga = SpacegroupAnalyzer(structure, symprec=symprec)
        structure = sga.get_primitive_standard_structure()

    st = CubicSupercellTransformation(
        min_atoms,
        max_atoms,
        num_nn_dists,
        force_diagonal_transformation=force_diagonal_transformation
    )
    supercell = st.apply_transformation(structure)

    perturbed_supercells, disp_dists, fws = get_perturbed_supercell_fws(
        supercell,
        common_settings,
        user_incar_settings=uis,
        min_displacement=min_displacement,
        max_displacement=max_displacement,
        num_displacements=num_displacements,
        min_random_distance=min_random_distance,
        n_cells_per_displacement_distance=n_cells_per_displacement_distance,
    )

    logger.debug("Using {} displacements".format(len(disp_dists)))
    logger.debug("Displacements: {}".format(disp_dists))

    # Collect force constants from the DB and output on cluster
    wf_meta = {"wf_uuid": str(uuid4()), "wf_name": "LatticeDynamicsWF"}
    csld_fw = Firework(
        CSLDForceConstantsToDB(
            db_file=common_settings["DB_FILE"],
            wf_uuid=wf_meta["wf_uuid"],
            parent_structure=structure,
            trans_mat=st.transformation_matrix.tolist(),
            supercell_structure=supercell,
            perturbed_supercells=perturbed_supercells,
            disps=disp_dists,
            force_diagonal_transformation=force_diagonal_transformation,
            first_pass=True,
            static_user_incar_settings=uis,
            env_vars=common_settings,
            dynamic_static_calcs=dynamic_static_calcs,
            do_shengbte=do_shengbte,
            shengbte_t_range=shengbte_t_range,
            shengbte_fworker=shengbte_fworker
        ),
        name="Compressed Sensing Lattice Dynamics",
        parents=fws[-len(perturbed_supercells):]
    )
    fws.append(csld_fw)

    formula = structure.composition.reduced_formula
    wf_name = "{} - compressed sensing lattice dynamics".format(formula)

    wf = Workflow(fws, name=wf_name)
    wf = add_additional_fields_to_taskdocs(
        wf,  {"wf_meta": wf_meta}, task_name_constraint="VaspToDb"
    )

    return wf


def get_perturbed_supercell_fws(
    supercell: Structure,
    common_settings: Optional[Dict] = None,
    user_incar_settings: Optional[Dict] = None,
    max_displacement: float = _MAX_DISPLACEMENT,
    min_displacement: float = _MIN_DISPLACEMENT,
    num_displacements: int = _NUM_DISPLACEMENTS,
    n_cells_per_displacement_distance: int = _CELLS_PER_DISPLACEMENT,
    min_random_distance: Optional[float] = _MIN_RANDOM_DISPLACEMENTS,
) -> List[Firework]:
    """
    Get static calculations to calculate the forces on each perturbed supercell.

    Args:
        supercell: Parent supercell structure.
        max_displacement: Maximum displacement distance for perturbing the
            sites in Angstroms.
        min_displacement: Minimum displacement distance for perturbing the
            sites in Angstroms.
        num_displacements: Number of unique displacement distances to
            try, uniformly distributed between ``min_displacement`` and
            ``max_displacement``.
        n_cells_per_displacement_distance: Number of perturbed
            structures to generate for each unique displacement distance.
        min_random_distance: If None (default), then for a given perturbed
            structure, all atoms will move the same distance from their original
            locations. If float, then for a given perturbed structure, the
            distances the atoms move will be uniformly distributed from a
            minimum distance of ``min_random_distance`` to one of the
            displacement distances uniformly sampled between
            ``min_displacement`` and ``max_displacement``.
        common_settings: Common settings dict
        user_incar_settings: User incar settings override.

    Returns:
        A list of static fireworks.
    """
    user_incar_settings = user_incar_settings if user_incar_settings else {}

    # Generate list of perturbed supercells
    perturbed_supercells, disp_dists = generate_perturbed_supercells(
        supercell,
        min_displacement=min_displacement,
        max_displacement=max_displacement,
        num_displacements=num_displacements,
        n_cells_per_displacement_distance=n_cells_per_displacement_distance,
        min_random_distance=min_random_distance
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
            name=name
        )
        static_fw.spec["displacement_value"] = dist
        fws.append(static_fw)

    return fws
