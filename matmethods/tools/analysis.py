# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import numpy as np

from pymatgen.analysis.eos import EOS

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'


def get_phonopy_gibbs(energies, volumes, force_constants, structure, t_min, t_step, t_max, mesh, eos):
    """
    Compute QHA gibbs free nergy using the phonopy interface.

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

    Returns:
        (numpy.ndarray, numpy.ndarray): Gibbs free energy , Temperature
    """
    try:
        from phonopy import Phonopy
        from phonopy.structure.atoms import Atoms as PhonopyAtoms
        from phonopy import PhonopyQHA
    except ImportError:
        import sys
        print("Install phonopy. Exiting.")
        sys.exit()

    phon_atoms = PhonopyAtoms(symbols=[str(s.specie) for s in structure],
                              scaled_positions=structure.frac_coords)
    phon_atoms.set_cell(structure.lattice.matrix)
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

    # quasi-harmonic approx
    phonopy_qha = PhonopyQHA(volumes, energies, eos=eos, temperatures=temperatures[0],
                             free_energy=np.array(free_energy).T, cv=np.array(cv).T,
                             entropy=np.array(entropy).T, t_max=np.max(temperatures[0]))

    # gibbs free energy and temperature
    max_t_index = phonopy_qha._qha._max_t_index
    G = phonopy_qha.get_gibbs_temperature()[:max_t_index]
    T = phonopy_qha._qha._temperatures[:max_t_index]
    return G, T


def get_debye_model_gibbs(energies, volumes, structure, t_min, t_step, t_max, eos):
    """
    Compute QHA gibbs free energy using the debye model. See the paper:
    http://doi.org/10.1016/j.comphy.2003.12.001

    Args:
        energies (list):
        volumes (list):
        structure (Structure):
        t_min (float): min temperature
        t_step (float): temperature step
        t_max (float): max temperature
        eos (str): equation of state used for fitting the energies and the volumes.
            options supported by pymatgen: "quadratic", "murnaghan", "birch", "birch_murnaghan",
            "pourier_tarantola", "vinet", "deltafactor"

    Returns:
        (numpy.ndarray, numpy.ndarray): Gibbs free energy , Temperature
    """
    temp = np.linspace(t_min, t_max, np.ceil((t_max-t_min)/t_step))
    mass = sum([e.atomic_mass for e in structure.species])
    natoms = structure.composition.num_atoms

    G = []
    T = []
    for t in temp:
        try:
            G_tmp = gibbs_minimizer(energies, volumes, mass, natoms, temperature=t, eos=eos)
        except:
            print("EOS fitting failed, so skipping this data point")
            continue
        G.append(G_tmp)
        T.append(t)
    return np.array(G), np.array(T)


def debye_integral(y, integrator):
    """
    Eq(5) in  doi.org/10.1016/j.comphy.2003.12.001

    Args:
        y (float)
        integrator: function from scipy.integrate

    Returns:
        float
    """
    if y < 155:
        return list(integrator(lambda x: x ** 3 / (np.exp(x) - 1), 0, y))[0] * 3. / y ** 3
    else:
        return 6.493939 * 3. / y ** 3  # limiting value: 6.49393940226684 from (wolfram alpha)


def A_vib(T, debye, natoms, integrator):
    """
    Eq(4) in  doi.org/10.1016/j.comphy.2003.12.001

    Args:
        T (float): temperature in K
        debye (float): debye temperature in K
        natoms (int): number of atoms
        integrator: function from scipy.integrate

    Returns:
        float: vibrational free energy
    """
    y = debye / T
    return 8.617332 * 1e-5 * natoms * T * (9. / 8. * y + 3 * np.log(1 - np.exp(-y))
                                           - debye_integral(y, integrator))


def debye_temperature_gibbs(volume, mass, natoms, B, poisson=0.25):
    """
    Shamelessly stole joey's implementation in pymatgen.

    Calculates the debye temperature accordings to the GIBBS formulation (in SI units)

    Args:
        volume (flaot)
        mass (float): total mass
        natoms (int): number of atoms
        B (float): bulk modulus
        poisson (flaot): poisson ratio

    Returns:
        debye temperature (in SI units)
     """
    avg_mass = 1.6605e-27 * mass / natoms
    t = poisson
    f = (3. * (2. * (2. / 3. * (1. + t) / (1. - 2. * t)) ** (1.5)
               + (1. / 3. * (1. + t) / (1. - t)) ** (1.5)) ** -1) ** (1. / 3.)
    return 2.9772e-11 * avg_mass ** (-1. / 2.) * (volume / natoms) ** (-1. / 6.) * f * B ** (0.5)


def gibbs_minimizer(energies, volumes, mass, natoms, temperature=298.0, pressure=0, poisson=0.25,
                    eos="murnaghan"):
    """

    Args:
        energies:
        volumes:
        mass (flaot): total mass
        natoms (int): number of atoms
        temperature (float): temperature in K
        pressure (flaot): pressure in GPa
        poisson (flaot): poisson ratio
        eos (str): name of the equation of state supported by pymatgen. See pymatgen.analysis.eos.py

    Returns:
        float: gibbs free energy at the given temperature and pressure minimized wrt volume.
    """
    try:
        from scipy.optimize import minimize
        from scipy.integrate import quadrature
    except ImportError:
        import sys
        print("Install scipy. Exiting.")
        sys.exit()

    integrator = quadrature

    eos = EOS(eos)
    eos_fit_1 = eos.fit(volumes, energies)
    G_V = []
    for i, v in enumerate(volumes):
        debye = debye_temperature_gibbs(v, mass, natoms, eos_fit_1.b0_GPa, poisson=poisson)
        G_V.append(energies[i] + pressure * v + A_vib(temperature, debye, natoms, integrator))
    eos_fit_2 = eos.fit(volumes, G_V)
    # get min vol
    params = eos_fit_2.eos_params.tolist()
    min_vol = minimize(eos_fit_2.func, min(volumes), tuple(params))
    return eos_fit_2.func(min_vol.x, *params)
