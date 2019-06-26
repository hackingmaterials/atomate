#!/usr/bin/env python3

import numpy as np

"""
This module defines a script that generates perturbed supercells from a supercell structure
"""

class CSLDPerturbedSupercellEnumerator:
    """
    Uses a supercell Structure, a transformation matrix, a min/max displacement,
    number of displacement values, and number of supercells per displacement
    to generate a list of randomly perturbed supercells from a given supercell.

        Args:
            original_supercell (Structure): supercell to displace, typically already relaxed
            max_displacement_val (float): maximum displacement value for perturbing the supercell (in Angstroms)
            min_displacement_val (float): minimum displacement value for perturbing the supercell (in Angstroms)
            num_displacements (int): number of unique displacement values to try, bounded by
                min_displacement and max_displacement. This argument is ignored if min_displacement
                is not passed as an argument.
            scs_per_displacement_val (int): number of perturbed supercells to generate for each unique
                displacement value.

            floor_displacement (Optional float): If None (default), then for a given perturbed supercell, all atoms will
                move the same distance from their original locations. If float, then for a given perturbed supercell,
                the distances that atoms move will be uniformly distributed from a minimum value of 'floor_displacement'
                to the displacement value 'displacement_val'.

    """

    def __init__(self, original_supercell,
                 max_displacement_val,
                 min_displacement_val,
                 num_displacements,
                 scs_per_displacement_val,
                 floor_displacement=None):

        self.original_supercell = original_supercell
        self.max_disp_val = max_displacement_val
        self.min_disp_val = min_displacement_val
        self.num_disps = num_displacements
        self.scs_per_disp_val = scs_per_displacement_val

        if floor_displacement is not None:
            self.floor_disp = float(floor_displacement)
        else:
            self.floor_disp = None

        self.disp_vals = np.linspace(min_displacement_val, max_displacement_val, num=num_displacements)
        self.perturbed_supercells = self.generate_transformations()


    def random_displacements(self, natom, rmax, rmin, dim=3):
        #*** Adapted from csld.util.mathtool

        # Generate matrix of size (natom, 3) where each row is a Gaussian random vector
        # If 'rmax'=None, then all vectors will have magnitude 'rmin'.
        # Else, magnitudes are uniformly distributed between 'rmin' and 'rmax'.

        dx = np.random.normal(size=(natom, dim)) #Gaussian sampled displacements. Matrix size: (natom, 3)
        dx_norms = np.linalg.norm(dx, axis=1)
        veclen = np.full(natom, rmax) if rmin is None else np.random.uniform(rmin, rmax, natom)

        if not 0 in dx_norms:
            return dx * (veclen / dx_norms)[:, None]
        else:
            self.random_displacements(natom, rmax, rmin, dim)


    def perturb_supercell(self, structure, max_disp, floor_disp):
        # *** Adapted from CSLD's 'polaron_main' file

        na = structure.num_sites
        max_disp = float(max_disp)
        if isinstance(max_disp, float):
            # from csld.util.mathtool import random_displacements
            dr = self.random_displacements(na, max_disp, floor_disp) #generate random displacements

            # Perturb structure with the random displacements
            for i in range(len(structure._sites)):
                structure.translate_sites([i], dr[i], frac_coords=False, to_unit_cell=True)
            return structure
        else:
            print('displacement entered is not a float.')


    def generate_transformations(self):
        perturbed_supercells = []
        for disp_val in self.disp_vals:
            for cell in range(self.scs_per_disp_val):
                sc = self.original_supercell.copy()
                sc = self.perturb_supercell(sc, disp_val, self.floor_disp)
                perturbed_supercells += [sc]
        return perturbed_supercells


####################################
# BELOW IS FOR DEBUGGING PURPOSES #
####################################
def main():
    # print('test')
    # natom = 10
    # dim = 3
    # rmin = 0
    # rmax = 0.1
    # dx = np.random.normal(size=(natom, dim))  # Gaussian sampled displacements. Matrix size: (natom, 3)
    # print(dx.shape)
    # veclen = np.full(natom, rmin) if rmax is None else np.random.uniform(rmin, rmax, natom)  # vector of length natom
    # dx *= (veclen / np.linalg.norm(dx, axis=1))[:, None]
    # test = (veclen / np.linalg.norm(dx, axis=1))[:, None]
    # print(test.shape)
    #
    # x = np.array([[1, 1, 1],
    #               [1, 1, 1],
    #               [1, 1, 1],
    #               [1, 1, 1]])
    # y = np.array([[1],
    #               [2],
    #               [3],
    #               [4]])
    # print(x*y)

    from atomate.vasp.workflows.base.csld_sc_gen import gen_scaling_matrix

    from pymatgen.ext.matproj import MPRester
    from pymatgen.analysis.structure_analyzer import get_max_bond_lengths
    from pymatgen.transformations.standard_transformations import SupercellTransformation

    mpr = MPRester(api_key='auNIrJ23VLXCqbpl')
    # structure = mpr.get_structure_by_material_id('mp-7814') #passed
    structure = mpr.get_structure_by_material_id('mp-1101039')  # passed
    # structure = mpr.get_structure_by_material_id('mp-20012')

    print('structure loaded')
    print(structure.num_sites)
    # print('structure')
    # print(structure)
    # structure.to(fmt='poscar', filename='POSCAR-tlbise2-scraped')

    # # Find Primitive
    # from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    # sga = SpacegroupAnalyzer(structure)
    # # print(sga.find_primitive().lattice)
    # sga.find_primitive().to(fmt='poscar', filename='POSCAR-tlbise2-scraped-prim')

    T = gen_scaling_matrix(structure, length_cutoff=5, min_atoms=150, max_atoms=1000)
    # ***my function, pymatgen's SupercellTransformation, and CSLD all yield different supercell lattice vectors

    print('~~~~~~~~~~~~~~~~~~~~~~~~~')
    print("FINISHED2")
    # print(T)

    # Pymatgen version
    superstructure = SupercellTransformation(T).apply_transformation(structure)

    enum = CSLDPerturbedSupercellEnumerator(superstructure,
                                            max_displacement_val=0.05,
                                            min_displacement_val=0.01,
                                            num_displacements=2,
                                            scs_per_displacement_val=1)
    print('perturbed supercells list:')
    print(enum.perturbed_supercells)


if __name__ == '__main__':
    main()