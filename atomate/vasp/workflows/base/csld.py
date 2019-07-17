# coding: utf-8

import os

import numpy as np
from uuid import uuid4

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.alchemy.materials import TransformedStructure
from pymatgen.transformations.advanced_transformations import (
    CubicSupercellTransformation,
    PerturbSitesTransformation
)
from fireworks import Workflow, Firework
from atomate.utils.utils import get_logger
from atomate.vasp.config import VASP_CMD, DB_FILE, ADD_WF_METADATA #                what is this?
from atomate.vasp.fireworks.core import StaticFW, OptimizeFW
from atomate.vasp.firetasks.parse_outputs import CSLDForceConstantsToDB
from atomate.vasp.powerups import (
    add_additional_fields_to_taskdocs
)

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
            symprec=0.1,
            min_atoms=-np.Inf,
            max_atoms=np.Inf,
            num_nn_dists=5,
            max_displacement=0.1,
            min_displacement=0.01,
            num_displacements=10,
            supercells_per_displacement_distance=1,
            min_random_distance=None,
            #csld input params here
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
            parent_structure (Structure):
            min_atoms (int):
            max_atoms (int):
            num_nn_dists (int or float):
            max_displacement (float)
        """
        self.uuid = str(uuid4())
        self.wf_meta = {
            "wf_uuid": self.uuid,
            "wf_name": self.__class__.__name__, #"CompressedSensingLatticeDynamicsWF"
        }

    # Create supercell
        sga = SpacegroupAnalyzer(parent_structure, symprec=symprec)
        self.parent_structure = sga.get_primitive_standard_structure()

        self.min_atoms = min_atoms
        self.max_atoms = max_atoms
        self.num_nn_dists = num_nn_dists
        supercell_transform = CubicSupercellTransformation(
            self.min_atoms,
            self.max_atoms,
            self.num_nn_dists,
        )
        # supercell (Structure)
        self.supercell = supercell_transform.apply_transformation(self.parent_structure)
        self.trans_mat = supercell_transform.trans_mat

    # Generate randomly perturbed supercells
        perturbed_supercells_transform = PerturbSitesTransformation(
            max_displacement,
            min_displacement,
            num_displacements,
            supercells_per_displacement_distance,
            min_random_distance
        )
        # list of perturbed supercell structures (list)
        self.perturbed_supercells = perturbed_supercells_transform.apply_transformation(self.supercell)
        # list of (non-unique) displacement values used in the perturbation (np.ndarray)
        self.disps = np.repeat(perturbed_supercells_transform.disps,
                               supercells_per_displacement_distance)

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

        #TODO: Move this relaxation section elsewhere
        # relax_user_incar_settings = {"EDIFF": 1e-8,
        #                              "EDIFFG": -1e-5,
        #                              }
        # relax_vis = MPRelaxSet(self.parent_structure,
        #                        user_incar_settings=relax_user_incar_settings)
        # # relax
        # fws.append(
        #     OptimizeFW(
        #         self.parent_structure,
        #         vasp_input_set=relax_vis,
        #         vasp_cmd=c["VASP_CMD"],
        #         db_file=c["DB_FILE"],
        #         max_force_threshold=0.05, #idk, should i change this?
        #         half_kpts_first_relax=False, #idk what this is
        #         name="{} - CSLD relax parent".format(self.parent_structure.composition.reduced_formula),
        #     )
        # )
        #####################

        static_user_incar_settings = {"ADDGRID": True,
                                      # Fast Fourier Transform grid
                                      "LCHARG": False,
                                      "ENCUT": 700,
                                      "EDIFF": 1e-7,  # may need to tune this
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
            fws.append(StaticFW(
                perturbed_supercell,
                vasp_input_set=static_vis,
                vasp_cmd=c["VASP_CMD"],
                db_file=c["DB_FILE"],
                name=name + " static"
            ))

        print('DISPS')
        print(self.disps)
        # Collect force constants from the DB and output on cluster
        csld_fw = Firework(
            CSLDForceConstantsToDB(
                db_file=c["DB_FILE"], # wot
                wf_uuid=self.uuid,
                name='CSLDForceConstantsToDB',
                parent_structure=self.parent_structure,
                trans_mat=self.trans_mat,
                supercell_structure=self.supercell,
                perturbed_supercells=self.perturbed_supercells,
                disps=self.disps
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
    from atomate.vasp.powerups import add_tags, set_execution_options

    #get a structure
    mpr = MPRester(api_key='auNIrJ23VLXCqbpl')
    structure = mpr.get_structure_by_material_id('mp-149')

    # prim_structure = structure.get_primitive_structure() #this method is bad

    # CONVENTIONAL PRIMITIVE CELL WITH SYMMETRIZATION
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    sga = SpacegroupAnalyzer(structure, symprec=0.1)
    prim = sga.get_primitive_standard_structure()
    prim.to("poscar", filename="POSCAR-pmg_prim_si")

    # CSLD'S POLARON MAIN-GENERATED PRIMITIVE CELL
    # from pymatgen.core.structure import Structure
    # prim = Structure.from_file('POSCAR_csld_primitivized')

    csld_class = CompressedSensingLatticeDynamicsWF(prim, max_displacement=0.05,
                                                    min_displacement=0.01,
                                                    num_displacements=5)
    print("uuid")
    print(csld_class.uuid)

    wf = csld_class.get_wf()
    print("trans mat")
    print(csld_class.trans_mat)
    csld_class.supercell.to("poscar", filename="POSCAR-csld_super_si")

    wf = add_tags(wf, ['csld', 'v1', 'rees', 'sga primitizer'])
    wf = set_execution_options(wf, fworker_name="rees_the_fire_worker") #need this to run the fireworks
    print(wf)

    # lpad = LaunchPad.auto_load()
    # lpad.add_wf(wf)

    # [[4.  0. - 2.]
    #  [-1.  4. - 1.]
    # [0.
    # 0.
    # 4.]]
