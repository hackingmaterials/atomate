#!/usr/bin/env python3

import math

import numpy as np

from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.structure_analyzer import get_max_bond_lengths
from pymatgen.transformations.standard_transformations import SupercellTransformation

"""
This module defines the Compressed Sensing Lattice Dynamics (CSLD) workflow.
"""

def gen_scaling_matrix(structure=None, max_atoms=np.Inf, min_atoms=-np.Inf, length_cutoff=5):
    """
    Returns the transformation matrix of a Pymatgen structure into a suitable supercell for CSLD

    Args:
        structure (Structure): input structure.
        max_atoms (int): maximum number of atoms allowed in the supercell
        min_atoms (int): minimum number of atoms allowed in the supercell
        length_cutoff (int): number of multiples of the unit cell nearest neighbor distance to
            force all directions of the supercell to be at least as large

    Returns:
        Transformation matrix (numpy array)
    """

    if not structure:
        print('No structure was passed into gen_scaling_matrix()')
    else:

        # -------- Store relevant structure data ----------
        latVecs = structure.lattice.matrix #lattice vectors
        bondLengths = get_max_bond_lengths(structure) #dictionary of bond lengths
        bondMatrix = structure.distance_matrix #NxN matrix of bond distances (diagonal is 0)
        print(bondMatrix)
        print(structure.sites)
        np.fill_diagonal(bondMatrix, np.Inf)
        nnDist = np.amin(bondMatrix) #nearest neighbor bond length
        print("nn Dist: " + str(nnDist))

        # --------- Transform lattice vectors -------------
        SCnotFound = True  # boolean indicator for whether a sufficiently large supercell has been created
        hard_threshold = nnDist * length_cutoff
        print("hard_threshold: " + str(hard_threshold))
        target_threshold = hard_threshold  # target_threshold is used as the desired CUBE side lengths of
                                            # the supercell (SC)
        # note: if it is not possible to get a cubic SC, the real SC dimensions will be less
        #      than the target
        while SCnotFound:
            print('---------------------------------------------------------')
            print('target threshold: ' + str(target_threshold))
            target_supercell = np.eye(3, 3) * target_threshold

            T = np.linalg.inv(latVecs) @ target_supercell
            # T = target_supercell @ np.linalg.inv(latVecs)
            T = np.around(T)
            print("proposed sc.txt:")
            print(T)

            print("\nSupercell dimensions:")
            scBases = latVecs @ T
            print(scBases)
            print("Supercell volume:")
            vol = np.cross(scBases[0],scBases[1]).T @ scBases[2]
            print(vol)

            # -------------------- Check how many nearest neighbors the supercell allows ---------------------------
            def proj(b, a):
                # Returns vector projection of b onto a
                return (b.T @ (a / np.linalg.norm(a))) * (a / np.linalg.norm(a))

            scBasesNorms = [np.linalg.norm(scBases[0]), np.linalg.norm(scBases[1]), np.linalg.norm(scBases[2])]
            maxNorm = max(scBasesNorms)
            maxIndex = scBasesNorms.index(maxNorm)
            idx = list(range(3))
            idx.remove(maxIndex)
            a = scBases[maxIndex]
            b = scBases[idx[0]]
            c = scBases[idx[1]]
            projb_a = proj(b, a)
            projc_a = proj(c, a)

            if np.linalg.norm(projb_a) > np.linalg.norm(projc_a):
                length = np.linalg.norm(a - projb_a)
            else:
                length = np.linalg.norm(a - projc_a)
            width = math.sqrt(np.linalg.norm(b) ** 2 - np.linalg.norm(projb_a) ** 2)
            ab_normal = np.cross(a, b)  # get normal direction from AB plane
            height = np.linalg.norm(proj(c, ab_normal))  # project c onto AB plane normal

            cubeSide = min([length, width, height])
            print('smallest dimension of proposed SC: ' + str(cubeSide))

            superstructure = SupercellTransformation(T).apply_transformation(structure)
            num_at = superstructure.num_sites
            if cubeSide >= hard_threshold and num_at >= min_atoms and num_at <= max_atoms:
                print(T)
                print("FINISHED")
                return T
            else:
                # Increase threshold until
                target_threshold += 1


def parallelipipedVol(a=[1, 0, 0], b=[0, 1, 0], c=[0, 0, 1]):
    vol = np.cross(a,b).T @ c
    return vol

def main():
    mpr = MPRester(api_key='auNIrJ23VLXCqbpl')
    structure = mpr.get_structure_by_material_id('mp-7814')
    print('structure loaded')

    T = gen_scaling_matrix(structure, length_cutoff=4, min_atoms=150)
    #***my function, pymatgen's SupercellTransformation, and CSLD all yield different supercell lattice vectors

    print('~~~~~~~~~~~~~~~~~~~~~~~~~')
    print("FINISHED2")
    # print(T)

    # Pymatgen version
    superstructure = SupercellTransformation(T).apply_transformation(structure)
    print(superstructure.lattice)
    print(superstructure.volume)
    print(parallelipipedVol(superstructure.lattice.matrix[0],
                            superstructure.lattice.matrix[1],
                            superstructure.lattice.matrix[2]))
    print(superstructure.num_sites)
    superstructure.to(fmt='poscar', filename='POSCAR-scTest')

    # Find Primitive
    # from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    # sga = SpacegroupAnalyzer(superstructure)
    # print(sga.find_primitive().lattice)

    # # CSLD version
    # # *** Can't call CSLD Structure member functions because structure scraped is a Pymatgen Structure ***
    # from csld.structure import Structure as CSLDStructure
    # from csld.lattice import Lattice as CSLDLattice
    # from csld.coord_utils import supercell_latticepoints
    #
    # def generate_supercell(CSLDStructure, scaling_matrix, scref=None):
    #     """
    #     Create a supercell.
    #
    #     Args:
    #         scaling_matrix: A scaling matrix for transforming the lattice
    #             vectors. Has to be all integers. Several options are possible:
    #
    #             a. A full 3x3 scaling matrix defining the linear combination
    #                the old lattice vectors. E.g., [[2,1,0],[0,3,0],[0,0,
    #                1]] generates a new structure with lattice vectors a' =
    #                2a + b, b' = 3b, c' = c where a, b, and c are the lattice
    #                vectors of the original structure.
    #             b. An sequence of three scaling factors. E.g., [2, 1, 1]
    #                specifies that the supercell should have dimensions 2a x b x
    #                c.
    #             c. A number, which simply scales all lattice vectors by the
    #                same factor.
    #     """
    #     scmat = np.array(scaling_matrix, np.int16)
    #     if scmat.shape != (3, 3):
    #         scmat= np.array(scmat* np.eye(3), np.int16)
    #     n_cell=int(round(np.linalg.det(scmat)))
    #     old_lattice = CSLDStructure._lattice
    #     new_lattice = CSLDLattice(np.dot(scmat, old_lattice.matrix))
    #     tvects = supercell_latticepoints(scmat)
    #     inv=np.linalg.inv(scmat)
    #     if scref is None:
    #         sc_ref= supercell_latticepoints(scmat)
    #     else:
    #         sc_ref= scref
    #     return CSLDStructure(CSLDLattice(np.dot(scmat, CSLDStructure._lattice.matrix)),
    #        [s.species_and_occu for s in CSLDStructure for _ in range(n_cell)],
    #        (CSLDStructure.frac_coords[:,None,:]+sc_ref[None,:,:]).reshape((-1,3)).dot(inv),
    #        coords_are_cartesian=False, to_unit_cell=True,
    #        site_properties_T=[s.properties for s in CSLDStructure for _ in range(n_cell)],
    #        intensive_properties=CSLDStructure.intensive_properties,extensive_properties=
    #                        {k:v*CSLDStructure.n_cell for k,v in CSLDStructure.extensive_properties.items()})
    #
    # # csld_structure = csld.structure.Structure(structure.lattice, structure.species, structure.lattice.get_fractional_coords())
    # csld_supercell = SupercellStructure.from_scmat(csld_structure, T)
    # print("here:")
    # print(csld_supercell)

if __name__ == '__main__':
    main()




