#!/usr/bin/env python3

import os
import sys

import numpy as np

"""
This module defines a script that generates perturbed supercells from a supercell structure
"""

class CSLDPerturbedSupercellEnumerator:
    """
    Uses a supercell Structure, a transformation matrix, a min/max displacement,
    number of displacement values, and number of supercells per displacement
    to generate a list of randomly perturbed supercells from a given supercell.
    """

    def __init__(self, original_supercell,
                 transformation_matrix,
                 min_displacement,
                 max_displacement,
                 num_displacements,
                 scs_per_displacement):

        self.supercell = original_supercell
        self.trans_mat = transformation_matrix
        self.min_disp = min_displacement
        self.max_disp = max_displacement
        self.num_disps = num_displacements
        self.scs_per_disp = scs_per_displacement
        self.perturbed_supercells = []


    def random_displacements(self, natom, rmax, rmin=None, dim=3):
        #*** Adapted from csld.util.mathtool

        # Generate matrix of size (natom, 3) where each row is a Gaussian random vector
        # If 'rmax'=None, then all vectors will have magnitude 'rmin'.
        # Else, magnitudes are uniformly distributed between 'rmin' and 'rmax'.

        dx = np.random.normal(size=(natom, dim)) #Gaussian sampled displacements. Matrix size: (natom, 3)
        dx_norms = np.linalg.norm(dx, axis=1)
        veclen = np.full(natom, rmax) if rmax is None else np.random.uniform(rmin, rmax, natom)

        if not 0 in dx_norms:
            return dx * (veclen / dx_norms)[:, None]
        else:
            self.random_displacements(natom, rmax, rmin, dim)


    def perturb_supercell(self, structure, disp):
        # *** Adapted from CSLD's 'polaron_main' file

        na = structure.num_sites
        disp = float(disp)
        if isinstance(disp, float):
            # from csld.util.mathtool import random_displacements
            dr = self.random_displacements(na, disp) #generate random displacements

            # Perturb structure with the random displacements
            for i in range(len(structure._sites)):
                structure.translate_sites([i], dr[i], frac_coords=False, to_unit_cell=True)
        else:
            print('displacement entered is not a float.')


    def generate_transformations(self):
        for disp_val in range(self.num_disps):
            for cell in range(self.scs_per_disp):
                sc = self.original_supercell.copy()
                self.perturb_supercell(sc, self.max_disp)
                self.perturbed_supercells += [sc]


####################################
# BELOW IS FOR DEBUGGING PURPOSES #
####################################
def main():
    print('test')
    natom = 10
    dim = 3
    rmin = 0
    rmax = 0.1
    dx = np.random.normal(size=(natom, dim))  # Gaussian sampled displacements. Matrix size: (natom, 3)
    print(dx.shape)
    veclen = np.full(natom, rmin) if rmax is None else np.random.uniform(rmin, rmax, natom)  # vector of length natom
    dx *= (veclen / np.linalg.norm(dx, axis=1))[:, None]
    test = (veclen / np.linalg.norm(dx, axis=1))[:, None]
    print(test.shape)

    x = np.array([[1, 1, 1],
                  [1, 1, 1],
                  [1, 1, 1],
                  [1, 1, 1]])
    y = np.array([[1],
                  [2],
                  [3],
                  [4]])
    print(x*y)



if __name__ == '__main__':
    main()