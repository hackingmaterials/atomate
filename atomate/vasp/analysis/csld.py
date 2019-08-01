# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import numpy as np

from pymatgen.transformations.standard_transformations import PerturbStructureTransformation

__author__ = 'Rees Chang'
__email__ = 'rc564@cornell.edu'


def generate_perturbed_supercells(supercell,
                                  min_displacement=0.01, max_displacement=0.1,
                                  num_displacements=10,
                                  supercells_per_displacement_distance=1,
                                  min_random_distance=None):
    """
    Generate a list of supercells with perturbed atomic sites.

    Args:
        supercell (Structure): original structure whose atomic sites are to be
            perturbed
        max_displacement (float): maximum displacement distance for
            perturbing the structure (Angstroms)
        min_displacement (float): minimum displacement distance for
            perturbing the structure (Angstroms)
        num_displacements (int): number of unique displacement distances to
            try, uniformly distributed between 'min_displacement' and
            'max_displacement'.
        structures_per_displacement_distance (int): number of perturbed
            structures to generate for each unique displacement distance.
        min_random_distance (Optional float): If None (default), then for a
            given perturbed structure, all atoms will move the same distance
            from their original locations. If float, then for a given
            perturbed structure, the distances that atoms move will be
            uniformly distributed from a minimum distance of
            'min_random_distance' to one of the displacement distances
            uniformly sampled between 'min_displacement' and
            'max_displacement'.
    Returns:
        (List of randomly displaced structures (List of Structures),
         List of corresponding displacements)
    """
    displacements = np.repeat(
        np.linspace(min_displacement, max_displacement, num_displacements),
        supercells_per_displacement_distance
    )
    # self.min_random_distance = min_random_distance
    perturbed_supercells = []  # list of perturbed supercell structures
    for displacement in displacements:
        perturb_structure_transformer = PerturbStructureTransformation(
            displacement,
            min_random_distance)
        perturbed_supercell = perturb_structure_transformer.apply_transformation(
            supercell)
        perturbed_supercells += [perturbed_supercell]
    return perturbed_supercells, displacements
