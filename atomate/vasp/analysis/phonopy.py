# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import numpy as np

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

# TODO: @matk86 - unit tests?

def get_phonopy_gibbs(energies, volumes, force_constants, structure, t_min, t_step, t_max, mesh,
                      eos, pressure=0):
    """
    Compute QHA gibbs free energy using the phonopy interface.

    Args:
        energies (list):
        volumes (list):
        force_constants (list):
        structure (Structure):
        t_min (float): min temperature
        t_step (float): temperature step
        t_max (float): max temperature
        mesh (list/tuple): reciprocal space density
        eos (str): equation of state used for fitting the energies and the volumes.
            options supported by phonopy: vinet, murnaghan, birch_murnaghan
        pressure (float): in GPa, optional.

    Returns:
        (numpy.ndarray, numpy.ndarray): Gibbs free energy, Temperature
    """

    # quasi-harmonic approx
    phonopy_qha = get_phonopy_qha(energies, volumes, force_constants, structure, t_min, t_step,
                                  t_max, mesh, eos, pressure=pressure)

    # gibbs free energy and temperature
    max_t_index = phonopy_qha._qha._max_t_index
    G = phonopy_qha.get_gibbs_temperature()[:max_t_index]
    T = phonopy_qha._qha._temperatures[:max_t_index]
    return G, T


def get_phonopy_qha(energies, volumes, force_constants, structure, t_min, t_step, t_max, mesh, eos,
                      pressure=0):
    """
    Return phonopy QHA interface.

    Args:
        energies (list):
        volumes (list):
        force_constants (list):
        structure (Structure):
        t_min (float): min temperature
        t_step (float): temperature step
        t_max (float): max temperature
        mesh (list/tuple): reciprocal space density
        eos (str): equation of state used for fitting the energies and the volumes.
            options supported by phonopy: vinet, murnaghan, birch_murnaghan
        pressure (float): in GPa, optional.

    Returns:
        PhonopyQHA
    """
    from phonopy import Phonopy
    from phonopy.structure.atoms import Atoms as PhonopyAtoms
    from phonopy import PhonopyQHA
    from phonopy.units import EVAngstromToGPa

    phon_atoms = PhonopyAtoms(symbols=[str(s.specie) for s in structure],
                              scaled_positions=structure.frac_coords,
                              cell=structure.lattice.matrix)
    scell = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    phonon = Phonopy(phon_atoms, scell)
    # compute the required phonon thermal properties
    temperatures = []
    free_energy = []
    entropy = []
    cv = []
    for f in force_constants:
        phonon.set_force_constants(-np.array(f))
        phonon.set_mesh(list(mesh))
        phonon.set_thermal_properties(t_step=t_step, t_min=t_min, t_max=t_max)
        t, g, e, c = phonon.get_thermal_properties()
        temperatures.append(t)
        free_energy.append(g)
        entropy.append(e)
        cv.append(c)

    # add pressure contribution
    energies = np.array(energies) + np.array(volumes) * pressure / EVAngstromToGPa
    # quasi-harmonic approx
    return PhonopyQHA(volumes, energies, eos=eos, temperatures=temperatures[0],
                      free_energy=np.array(free_energy).T, cv=np.array(cv).T,
                      entropy=np.array(entropy).T, t_max=np.max(temperatures[0]))


def get_phonopy_thermal_expansion(energies, volumes, force_constants, structure, t_min, t_step,
                                  t_max, mesh, eos, pressure=0):
    """
    Compute QHA thermal expansion coefficient using the phonopy interface.

    Args:
        energies (list):
        volumes (list):
        force_constants (list):
        structure (Structure):
        t_min (float): min temperature
        t_step (float): temperature step
        t_max (float): max temperature
        mesh (list/tuple): reciprocal space density
        eos (str): equation of state used for fitting the energies and the volumes.
            options supported by phonopy: vinet, murnaghan, birch_murnaghan
        pressure (float): in GPa, optional.

    Returns:
        (numpy.ndarray, numpy.ndarray): thermal expansion coefficient, Temperature
    """

    # quasi-harmonic approx
    phonopy_qha = get_phonopy_qha(energies, volumes, force_constants, structure, t_min, t_step,
                                  t_max, mesh, eos, pressure=pressure)

    # thermal expansion coefficient and temperature
    max_t_index = phonopy_qha._qha._max_t_index
    alpha = phonopy_qha.get_thermal_expansion()[:max_t_index]
    T = phonopy_qha._qha._temperatures[:max_t_index]
    return alpha, T
