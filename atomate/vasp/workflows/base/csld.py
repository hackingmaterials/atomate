# coding: utf-8

import os

import numpy as np
from uuid import uuid4

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.alchemy.materials import TransformedStructure
from pymatgen.transformations.standard_transformations import \
    PerturbStructureTransformation
from pymatgen.transformations.advanced_transformations import (
    CubicSupercellTransformation
)
from fireworks import Workflow, Firework
from atomate.utils.utils import get_logger
from atomate.vasp.config import VASP_CMD, DB_FILE, ADD_WF_METADATA #                what is this?
from atomate.vasp.fireworks.core import StaticFW, OptimizeFW
from atomate.vasp.firetasks.parse_outputs import CSLDForceConstantsToDB
from atomate.vasp.powerups import (
    add_additional_fields_to_taskdocs
)
from atomate.vasp.analysis.csld import generate_perturbed_supercells

from pymatgen.io.vasp.sets import MPStaticSet, MPRelaxSet

# NEW THINGS TO INSTALL?
# from configparser import ConfigParser
# from scripts import csld_main_rees #csld>scripts

logger = get_logger(__name__)

__author__ = "Rees Chang"
__email__ = "rc564@cornell.edu"
__date__ = "July 2019"

__csld_wf_version__ = 1.0

class CompressedSensingLatticeDynamicsWF:
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
            shengbte_t_range=False,
            shengbte_fworker=None
    ):
        """
        This workflow will use compressed sensing lattice dynamics (CSLD)
        (doi: 10.1103/PhysRevLett.113.185501) to generate interatomic force
        constants from an input structure and output a summary to a database.

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
            shengbte_t_range (bool): If False (default), calculate the lattice
                thermal conductivity with ShengBTE at 300 K. If True, do the
                calculation from 100-1000 K with 100 K steps (much slower).
            shengbte_fworker (None or str): If None, the ShengBTE firework's
                fworker will be set to all the previous fireworks' fworker. If
                str, the ShengBTE firework's fworker will be set to
                shengbte_fworker.
        """
        # if force_diagonal_transformation is True and num_nn_dists==5:
        #     num_nn_dists = 6

        self.uuid = str(uuid4())
        self.wf_meta = {
            "wf_uuid": self.uuid,
            "wf_name": self.__class__.__name__, #"CompressedSensingLatticeDynamicsWF"
        }

    # Create supercell
        self.parent_structure = parent_structure
        if symmetrize:
            sga = SpacegroupAnalyzer(parent_structure, symprec=symprec)
            self.parent_structure = sga.get_primitive_standard_structure()

        self.min_atoms = min_atoms
        self.max_atoms = max_atoms
        self.num_nn_dists = num_nn_dists
        supercell_transform = CubicSupercellTransformation(
            self.min_atoms,
            self.max_atoms,
            self.num_nn_dists,
            force_diagonal_transformation=force_diagonal_transformation
        )
        # supercell (Structure)
        self.supercell = supercell_transform.apply_transformation(self.parent_structure)
        self.trans_mat = supercell_transform.trans_mat.tolist()
        self.supercell_smallest_dim = supercell_transform.smallest_dim

    # Generate list of perturbed supercells
        # self.perturbed_supercells: list of perturbed supercells (list of Structures)
        # self.disps: list of (non-unique) displacement values used in the perturbation (np.array)
        self.perturbed_supercells, self.disps = generate_perturbed_supercells(
            self.supercell,
            min_displacement=min_displacement,
            max_displacement=max_displacement,
            num_displacements=num_displacements,
            supercells_per_displacement_distance=supercells_per_displacement_distance,
            min_random_distance=min_random_distance
            )

        self.shengbte_t_range = shengbte_t_range
        self.shengbte_fworker = shengbte_fworker

    def get_wf(
            self,
            c=None
    ):
        fws = []
        c = c or {"VASP_CMD": VASP_CMD, "DB_FILE": DB_FILE}

        def _add_metadata(structure):
            """
            Add metadata for easy querying from the database later.
            """
            return TransformedStructure(
                structure, other_parameters={"wf_meta": self.wf_meta}
            )


        static_user_incar_settings = {"ADDGRID": True,
                                      # Fast Fourier Transform grid
                                      "LCHARG": False,
                                      "ENCUT": 700,
                                      "EDIFF": 1e-8,  # may need to tune this
                                      "PREC": 'Accurate',
                                      "LAECHG": False,
                                      "LREAL": False,
                                      "LASPH": True}
        static_user_incar_settings.update(c.get("user_incar_settings", {}))

        for idx, perturbed_supercell in enumerate(self.perturbed_supercells):
            # Run static calculations on the perturbed supercells to compute forces on each atom
            name = "perturbed supercell, idx: {}, disp_val: {:.3f},".format(idx, self.disps[idx])

            static_vis = MPStaticSet(perturbed_supercell,
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

        print('DISPS')
        print(self.disps)
        # Collect force constants from the DB and output on cluster
        csld_fw = Firework(
            CSLDForceConstantsToDB(
                db_file=c["DB_FILE"], # wot
                wf_uuid=self.uuid,
                parent_structure=self.parent_structure,
                trans_mat=self.trans_mat,
                supercell_structure=self.supercell,
                supercell_smallest_dim=self.supercell_smallest_dim,
                perturbed_supercells=self.perturbed_supercells,
                disps=self.disps,
                first_pass=True,
                static_user_incar_settings=static_user_incar_settings,
                env_vars=c,
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

        wf = add_additional_fields_to_taskdocs(wf,
                                               {"wf_meta": self.wf_meta},
                                               task_name_constraint="VaspToDb"
                                               #may need to change this to "CSLDForceConstantsToDB"?
                                               )
        #tag =   #insert anything relevant to every firework in the workflow
        # wf = add_tags(wf, [tag, <insert whatever string you'd like>])
        return wf


# SCRIPT FOR CREATING THE WORKFLOW AND ADDING IT TO THE DATABASE
if __name__ == "__main__":

    from fireworks import LaunchPad
    from pymatgen.ext.matproj import MPRester
    from atomate.vasp.powerups import add_tags, set_execution_options, add_modify_incar

    #get a structure
    # mpr = MPRester(api_key='auNIrJ23VLXCqbpl')
    # # structure = mpr.get_structure_by_material_id('mp-149') #Si
    # structure = mpr.get_structure_by_material_id('mp-1101039')
    # test_trans_mat = np.linalg.inv(structure.lattice.matrix) @ (np.eye(3,3) * 20)
    # print(test_trans_mat)
    # print(np.linalg.det(test_trans_mat))
    # structure_test = structure * test_trans_mat
    # structure_test.to("poscar", filename="POSCAR-noninteger_trans_test")

    # # CONVENTIONAL PRIMITIVE CELL WITH SYMMETRIZATION
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    # sga = SpacegroupAnalyzer(structure, symprec=0.1)
    # prim = sga.get_primitive_standard_structure()

    # CSLD'S POLARON MAIN-GENERATED PRIMITIVE CELL
    from pymatgen.core.structure import Structure
    # prim = Structure.from_file('POSCAR_csld_primitivized')

    prim = Structure.from_file('POSCAR-Si555')
    # prim = Structure.from_file('POSCAR-Mg3Sb2_save')
    # print(prim)

    csld_class = CompressedSensingLatticeDynamicsWF(prim,
                                                    symmetrize=False,
                                                    num_nn_dists=6,
                                                    num_displacements=10,
                                                    supercells_per_displacement_distance=1,
                                                    force_diagonal_transformation=True
                                                    )
    print("uuid")
    print(csld_class.uuid)

    wf = csld_class.get_wf()
    print("trans mat")
    print(csld_class.trans_mat)
    print(type(csld_class.trans_mat))
    print(csld_class.supercell_smallest_dim)
    print(csld_class.supercell.num_sites)
    csld_class.supercell.to("poscar", filename="SPOSCAR-csld_super_Si")
    # csld_class.supercell.to("poscar", filename="SPOSCAR-Mg3Sb2")

    wf = add_tags(wf, ['csld', 'v1', 'rees',
                       'pre-relaxed si', 'diagonal supercell',
                       'not symmetrized', '555'])
    wf = set_execution_options(wf, fworker_name="rees_the_fire_worker") #need this to run the fireworks
    wf = add_modify_incar(wf,
                          modify_incar_params={'incar_update': {'ENCUT': 500,
                                                                'ISPIN': 1}})
    # print(wf)

    # lpad = LaunchPad.auto_load()
    # lpad.add_wf(wf)

# uuid
# 0ff0d430-b43f-4cb5-ac8b-55c465b7867c
# DISPS
# [0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 ]
# trans mat
# [[-3.  3.  3.]
#  [ 3. -3.  3.]
#  [ 3.  3. -3.]]
# 23.745154105620003
# 432

#################custom csld_main below##################
def csld_main(options, settings):
    """
    Runs CSLD minimization.

    Changes from original version:
        - Made 'prim' an argument in 'phonon_step()'
        - Moved execution files to this main() function to be called from
          atomate
        - Rewrote 'add_common_parameter' in 'common_main' to treat 'options' as
          a dictionary instead of ArgumentParser
    """
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib.backends.backend_pdf import PdfPages
    import atexit
    import csld
    from csld.symmetry_structure import SymmetrizedStructure

    from csld.lattice_dynamics import init_ld_model

    from csld.common_main import upon_exit, \
        init_training
    from csld.csld_main_functions import phonon_step, \
        save_pot, predict, fit_data


    freq_matrix = None #Rees
    pdfout = PdfPages(
        options['pdfout'].strip()) if options['pdfout'].strip() else None
    atexit.register(upon_exit, pdfout)

    prim = SymmetrizedStructure.init_structure(settings['structure'],
                                               options['symm_step'],
                                               options['symm_prim'],
                                               options['log_level'])
    model = init_ld_model(prim, settings['model'], settings[
        'LDFF'] if 'LDFF' in settings.sections() else {}, options['clus_step'],
                          options['symC_step'], options['ldff_step'])
    Amat, fval = init_training(model, settings['training'], options['train_step'])
    ibest, solutions, rel_err = fit_data(model, Amat, fval, settings['fitting'],
                                options['fit_step'], pdfout)
    if settings.has_section('phonon'):
        phonon, freq_matrix = phonon_step(model, solutions, settings['phonon'],
                             options['phonon_step'], pdfout, prim, return_eigen=True)
    if settings.has_section('export_potential'):
        save_pot(model, solutions[ibest], settings['export_potential'],
                 options['save_pot_step'], phonon)
    if settings.has_section('prediction'):
        predict(model, solutions, settings['prediction'], options['pred_step'])

    #OUTPUT
    # freq_matrix is (nbands, nkpoints) = frequencies. Check for negative entries
    # rel_err is cross validation error in percent
    return rel_err, freq_matrix #also want to return force constants



