# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import numpy as np

from pymatgen.analysis.eos import EOS

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'


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


def get_debye_model_gibbs(energies, volumes, structure, t_min, t_step, t_max, eos, pressure=0):
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
        pressure (float): in GPa, optional.

    Returns:
        (numpy.ndarray, numpy.ndarray): Gibbs free energy , Temperature
        Note: The data points for which the equation of state fitting fails are skipped.
    """
    temp = np.linspace(t_min, t_max, np.ceil((t_max-t_min)/t_step)+1)
    mass = sum([e.atomic_mass for e in structure.species])
    natoms = structure.composition.num_atoms

    G = []
    T = []
    for t in temp:
        try:
            G_tmp = gibbs_minimizer(energies, volumes, mass, natoms, temperature=t, eos=eos,
                                    pressure=pressure)
        except:
            print("EOS fitting failed, so skipping this data point")
            continue
        G.append(G_tmp)
        T.append(t)
    return G, T


def debye_integral(y, integrator):
    """
    Debye integral. Eq(5) in  doi.org/10.1016/j.comphy.2003.12.001

    Args:
        y (float)
        integrator: function from scipy.integrate

    Returns:
        float
    """
    # floating point limit is reached around y=155, so values beyond that are set to the
    # limiting value(T-->0, y --> \infty) of 6.4939394 (from wolfram alpha).
    factor = 3. / y ** 3
    if y < 155:
        return list(integrator(lambda x: x ** 3 / (np.exp(x) - 1.), 0, y))[0] * factor
    else:
        return 6.493939 * factor


def A_vib(T, debye, natoms, integrator):
    """
    Vibrational free energy. Eq(4) in doi.org/10.1016/j.comphy.2003.12.001

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
    Calculates the debye temperature. Eq(6) in doi.org/10.1016/j.comphy.2003.12.001
    Thanks to Joey.

    Args:
        volume (float)
        mass (float): total mass
        natoms (int): number of atoms
        B (float): bulk modulus
        poisson (float): poisson ratio

    Returns:
        debye temperature (in SI units)
     """
    avg_mass = 1.6605e-27 * mass / natoms
    term1 = (2. / 3. * (1. + poisson) / (1. - 2. * poisson)) ** (1.5)
    term2 = (1. / 3. * (1. + poisson) / (1. - poisson)) ** (1.5)
    f = (3. / (2. * term1 + term2)) ** (1. / 3.)
    return 2.9772e-11 * (volume / natoms) ** (-1. / 6.) * f * np.sqrt(B/avg_mass)


def gibbs_minimizer(energies, volumes, mass, natoms, temperature=298.0, pressure=0, poisson=0.25,
                    eos="murnaghan"):
    """
    Fit the input energies and volumes to the equation of state to obtain the bulk modulus which is
    subsequently used to obtain the debye temperature. The debye temperature is then used to compute
    the  vibrational free energy and the gibbs free energy as a function of volume, temperature and
    pressure. A second fit is preformed to get the functional form of gibbs free energy:(G, V, T, P).
    Finally G(V, P, T) is minimized with respect to V and the optimum value of G evaluated at V_opt,
    G_opt(V_opt, T, P), is returned.

    Args:
        energies (list): list of energies
        volumes (list): list of volumes
        mass (float): total mass
        natoms (int): number of atoms
        temperature (float): temperature in K
        pressure (float): pressure in GPa
        poisson (float): poisson ratio
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

    # G(V, T, P)
    eos_fit_2 = eos.fit(volumes, G_V)
    params = eos_fit_2.eos_params.tolist()
    # G_opt(V_opt, T, P)
    return params[0]


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
    try:
        from phonopy import Phonopy
        from phonopy.structure.atoms import Atoms as PhonopyAtoms
        from phonopy import PhonopyQHA
        from phonopy.units import EVAngstromToGPa
    except ImportError:
        import sys
        print("Install phonopy. Exiting.")
        sys.exit()

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
