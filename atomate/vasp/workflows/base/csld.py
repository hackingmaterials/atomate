# coding: utf-8

import numpy as np
from uuid import uuid4

from pymatgen.transformations.advanced_transformations import (
    CubicSupercellTransformation)
from pymatgen.io.vasp.sets import MPStaticSet

from fireworks import Workflow, Firework
from atomate.utils.utils import get_logger
from atomate.vasp.config import VASP_CMD, DB_FILE
from atomate.vasp.fireworks.core import StaticFW
from atomate.vasp.firetasks.parse_outputs import CSLDForceConstantsToDB
from atomate.vasp.powerups import add_additional_fields_to_taskdocs
from atomate.vasp.analysis.csld import generate_perturbed_supercells

__author__ = "Rees Chang"
__email__ = "rc564@cornell.edu"
__date__ = "July 2019"
__csld_wf_version__ = 1.0

logger = get_logger(__name__)


class CompressedSensingLatticeDynamicsWF:
    """
    This workflow will use compressed sensing lattice dynamics (CSLD)
    (doi: 10.1103/PhysRevLett.113.185501) to generate interatomic force
    constants from an input structure and output a summary to a database.

    TODO: Implement 'tight_relaxation', 'dfpt_check', and 'nscf_check' options

    A summary of the workflow is as follows:
       1. Transform the input structure into a supercell
       2. Transform the supercell into a list of supercells with all atoms
           randomly perturbed from their original sites
       3. Run static VASP calculations on each perturbed supercell to
           calculate atomic forces.
       4. Aggregate the forces and conduct the CSLD minimization algorithm
           to compute interatomic force constants.
       5. Output the interatomic force constants to the database.

    Args:
       parent_structure (Structure): Structure to conduct CSLD on
       symmetrize (bool): If True, the parent_structure will be primitized
           before running CSLD. Else, the parent_structure will be used
           for CSLD as is.
       tight_relaxation (bool): If True, run a tight relaxation step on the
           parent_structure. Else, do not. TODO: IMPLEMENT
       dfpt_check (bool): If True, run a DFPT calculation at the gamma
           point to check for dynamic instabilities. TODO: IMPLEMENT
       nscf_check (bool): If True, run a nscf electronic bands calculation to
           determine if parent_structure is a metal or not. If metal, set
           ISMEAR=1 for the static perturbed supercell calculations. Else, set
           ISMEAR=0. TODO: IMPLEMENT
       symprec (float): symmetry precision parameter for
           SpacegroupAnalyzer's 'get_primitive_standard_structure()'
           function
       min_atoms (int): Minimum number of atoms to constrain the supercell
       max_atoms (int): Maximum number of atoms to constrain the supercell
       num_nn_dists (int or float): Number of nearest neighbor distances
           that the shortest direction of the supercell must be as long as
       force_diagonal_transformation (bool): If True, the supercell
           transformation will be constrained to have a diagonal
           transformation matrix. If False, the supercell transformation
           will not be constrained to be diagonal (resulting in a more
           cubic supercell).
       max_displacement (float): Maximum displacement distance for
           perturbing the supercell structure (Angstroms)
       min_displacement (float): Minimum displacement distance for
           perturbing the supercell structure (Angstroms)
       num_displacements (int): Number of unique displacement distances to
           generate uniformly between 'min_displacement' and
           'max_displacement'
       min_random_distance (Optional float): If None (default), then all
           atoms in the supercell will move the same distance from their
           original locations. If float, then for a given supercell, the
           distances that atoms move will be uniformly distributed from a
           minimum distance of 'min_random_distance' to one of the
           displacement distances uniformly sampled between
           'min_displacement' and 'max_displacement'.
       supercells_per_displacement_distance (int): number of supercells to
           generate for each unique displacement distance.
       dynamic_static_calcs (bool): If True and CSLD fails to return all real
            harmonic phonon frequencies on the first try, then the CSLD firework
            will dynamically create more static perturbed supercell calculations
            with larger displacements with the aim of handling materials with
            loosely bound atoms (e.g. rattling modes).
       do_shengbte (bool): If True and CSLD successfully returns all real
            harmonic phonon frequencies, then the CSLD firework will dynamically
            create a ShengBTE firework for calculating the lattice thermal
            conductivity matrix and storing it to a database.
       shengbte_t_range (bool): If False (default), calculate the lattice
           thermal conductivity with ShengBTE at 300 K. If True, do the
           calculation from 100-1000 K with 100 K steps (much slower).
       shengbte_fworker (None or str): If None, the ShengBTE firework's
           fworker will be set to all the previous fireworks' fworker. If
           str, the ShengBTE firework's fworker will be set to
           shengbte_fworker.
   """

    def __init__(
            self,
            parent_structure,
            symmetrize=False,
            symprec=0.1,
            min_atoms=-np.Inf,
            max_atoms=np.Inf,
            num_nn_dists=5,
            force_diagonal_transformation=True,
            max_displacement=0.1,
            min_displacement=0.01,
            num_displacements=10,
            min_random_distance=None,
            supercells_per_displacement_distance=1,
            dynamic_static_calcs=False,
            do_shengbte=False,
            shengbte_t_range=False,
            shengbte_fworker=None
    ):

        if force_diagonal_transformation is True and num_nn_dists == 5:
            num_nn_dists = 6

        # default name is CompressedSensingLatticeDynamicsWF
        self.uuid = str(uuid4())
        self.wf_meta = {
            "wf_uuid": self.uuid,
            "wf_name": self.__class__.__name__,
        }

        # Create supercell
        self.parent_structure = parent_structure
        if symmetrize:
            sga = SpacegroupAnalyzer(parent_structure, symprec=symprec)
            self.parent_structure = sga.get_primitive_standard_structure()

        self.force_diagonal_transformation = force_diagonal_transformation
        self.min_atoms = min_atoms
        self.max_atoms = max_atoms
        self.num_nn_dists = num_nn_dists
        supercell_transform = CubicSupercellTransformation(
            self.min_atoms,
            self.max_atoms,
            self.num_nn_dists,
            force_diagonal_transformation=self.force_diagonal_transformation
        )

        self.supercell = supercell_transform.apply_transformation(
            self.parent_structure)
        self.nn_dist = supercell_transform.nn_dist
        self.trans_mat = supercell_transform.trans_mat.tolist()
        self.supercell_smallest_dim = supercell_transform.smallest_dim

        # Generate list of perturbed supercells
        # perturbed_supercells: list of perturbed supercells (list[Structures])
        # disps: list of (non-unique) displacement values used in the
        # perturbation (np.array)
        self.perturbed_supercells, self.disps = generate_perturbed_supercells(
            self.supercell,
            min_displacement=min_displacement,
            max_displacement=max_displacement,
            num_displacements=num_displacements,
            supercells_per_displacement_distance=supercells_per_displacement_distance,
            min_random_distance=min_random_distance
            )

        self.dynamic_static_calcs = dynamic_static_calcs
        self.do_shengbte = do_shengbte
        self.shengbte_t_range = shengbte_t_range
        self.shengbte_fworker = shengbte_fworker

    def get_wf(self, c=None):
        fws = []
        c = c or {"VASP_CMD": VASP_CMD, "DB_FILE": DB_FILE}

        static_user_incar_settings = {
            "ADDGRID": True,  # Fast Fourier Transform grid
            "LCHARG": False,
            "ENCUT": 700,
            "EDIFF": 1e-8,  # may need to tune this
            "PREC": 'Accurate',
            "LAECHG": False,
            "LREAL": False,
            "LASPH": True}
        static_user_incar_settings.update(c.get("user_incar_settings", {}))

        for idx, perturbed_supercell in enumerate(self.perturbed_supercells):
            # Run static calculations on the perturbed supercells to compute
            # forces on each atom
            name = "perturbed supercell, idx: {}, disp_val: {:.3f},".format(
                idx, self.disps[idx])

            static_vis = MPStaticSet(
                perturbed_supercell,
                user_incar_settings=static_user_incar_settings)
            static_fw = StaticFW(
                perturbed_supercell,
                vasp_input_set=static_vis,
                vasp_cmd=c["VASP_CMD"],
                db_file=c["DB_FILE"],
                name=name + " static"
            )
            static_fw.spec["displacement_value"] = self.disps[idx]
            fws.append(static_fw)

        logger.debug("using {} displacements".format(len(self.disps)))
        logger.debug("displacements: {}".format(self.disps))

        # Collect force constants from the DB and output on cluster
        csld_fw = Firework(
            CSLDForceConstantsToDB(
                db_file=c["DB_FILE"],
                wf_uuid=self.uuid,
                parent_structure=self.parent_structure,
                trans_mat=self.trans_mat,
                supercell_structure=self.supercell,
                supercell_smallest_dim=self.supercell_smallest_dim,
                perturbed_supercells=self.perturbed_supercells,
                disps=self.disps,
                force_diagonal_transformation=self.force_diagonal_transformation,
                first_pass=True,
                static_user_incar_settings=static_user_incar_settings,
                env_vars=c,
                dynamic_static_calcs=self.dynamic_static_calcs,
                do_shengbte=self.do_shengbte,
                shengbte_t_range=self.shengbte_t_range,
                shengbte_fworker=self.shengbte_fworker
            ),
            name="Compressed Sensing Lattice Dynamics",
            parents=fws[-len(self.perturbed_supercells):]
        )
        fws.append(csld_fw)

        formula = self.parent_structure.composition.reduced_formula
        wf_name = "{} - compressed sensing lattice dynamics".format(formula)
        wf = Workflow(fws, name=wf_name)

        wf = add_additional_fields_to_taskdocs(
            wf,  {"wf_meta": self.wf_meta}, task_name_constraint="VaspToDb")

        return wf


# SCRIPT FOR CREATING THE WORKFLOW AND ADDING IT TO THE DATABASE
if __name__ == "__main__":

    from fireworks import LaunchPad
    from atomate.vasp.powerups import add_tags, set_execution_options, \
        add_modify_incar
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.core.structure import Structure

    prim = Structure.from_file('POSCAR-well_relaxed_KSnBi')

    csld_class = CompressedSensingLatticeDynamicsWF(
        prim,
        symmetrize=False,
        num_nn_dists=6,
        num_displacements=10,
        supercells_per_displacement_distance=1,
        force_diagonal_transformation=True,
        do_shengbte=True,
        shengbte_fworker="rees_the_fire_worker_haswell"
        )
    print("uuid")
    print(csld_class.uuid)

    wf = csld_class.get_wf()
    print("trans mat")
    print(csld_class.trans_mat)
    print("nn dist")
    print(csld_class.nn_dist)
    print("supercell shortest direction (Angstroms)")
    print(csld_class.supercell_smallest_dim)
    print("supercell number of atoms")
    print(csld_class.supercell.num_sites)
    csld_class.supercell.to("poscar", filename="SPOSCAR-KSnBi_diagonal")

    wf = add_tags(wf, ['csld', 'v1', 'rees',
                       'pre-relaxed ksnbi', 'diagonal supercell',
                       'not symmetrized', 'ismear manually set to 0'])
    wf = set_execution_options(wf, fworker_name="rees_the_fire_worker")
    wf = add_modify_incar(wf,
                          modify_incar_params={
                              'incar_update': {
                                  'ENCUT': 500,
                                  'ISMEAR': 0,
                                  'ISPIN': 1}})

    lpad = LaunchPad.auto_load()
    lpad.add_wf(wf)
