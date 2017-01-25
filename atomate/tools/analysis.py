# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

from collections import  defaultdict
import numpy as np
from scipy.integrate import quadrature

from pymatgen.analysis.eos import EOS

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'


class DebyeModelQHA(object):
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
    """
    def __init__(self, energies, volumes, structure, t_min, t_step, t_max, eos, pressure=0,
                 poisson=0.25):
        self.energies = energies
        self.volumes = volumes
        self.structure = structure
        self.temperature_min = t_min
        self.temperature_max = t_max
        self.temperature_step = t_step
        self.eos = eos
        self.pressure = pressure
        self.poisson = poisson
        self.mass = sum([e.atomic_mass for e in self.structure.species])
        self.natoms = self.structure.composition.num_atoms
        self.gibbs_free_energy = []
        self.temperatures = []
        self.optimum_volumes = []
        # fit E and V, get the bulk modulus
        eos = EOS(eos)
        eos_fit = eos.fit(volumes, energies)
        self.bulk_modulus = eos_fit.b0_GPa  # in GPa

    def compute_gibbs_free_energy(self):
        """
        Evaluate the gibbs free energy as a function of V,T and P i.e G(V, T, P) and
        minimize G(V, T, P) wrt V.

        Note: The data points for which the equation of state fitting fails are skipped.
        """
        temp = np.linspace(self.temperature_min, self.temperature_max,
                           np.ceil((self.temperature_max-self.temperature_min)/self.temperature_step)+1)
        for t in temp:
            try:
                G_opt, V_opt = self.minimizer(t)
            except:
                print("EOS fitting failed, so skipping this data point")
                continue
            self.gibbs_free_energy.append(G_opt)
            self.temperatures.append(t)
            self.optimum_volumes.append(V_opt)

    def minimizer(self, temperature):
        """
        1. Use the debye temperature to compute the  vibrational free energy.
        2. Compute the gibbs free energy as a function of volume, temperature and pressure.
        3. A second fit is preformed to get the functional form of gibbs free energy:(G, V, T, P).
        4. Finally G(V, P, T) is minimized with respect to V and the optimum value of G evaluated
            at V_opt, G_opt(V_opt, T, P) and V_opt are returned.

        Args:
            temperature (float): temperature in K

        Returns:
            float, float: gibbs free energy at the given temperature and pressure minimized wrt
                volume.
        """
        G_V = []  # G for each volume

        # G = E + PV + A_vib
        for i, v in enumerate(self.volumes):
            G_V.append(self.energies[i] +
                       self.pressure * v +
                       self.vibrational_free_energy(temperature, v))

        # fit G(V, T, P)
        eos_fit_2 = self.eos.fit(self.volumes, G_V)
        params = eos_fit_2.eos_params.tolist()
        # G_opt(V_opt, T, P), V_opt
        return params[0], params[-1]

    def vibrational_free_energy(self, temperature, volume):
        """
        Vibrational Helmholtz free energy, A_vib(V, T).
        Eq(4) in doi.org/10.1016/j.comphy.2003.12.001

        Args:
            temperature (float): temperature in K
            volume (float)

        Returns:
            float: vibrational free energy
        """
        y = self.debye_temperature(volume) / temperature
        return 8.617332 * 1e-5 * self.natoms * temperature * (9. / 8. * y + 3 * np.log(1 - np.exp(-y))
                                                              - self.debye_integral(y))

    def vibrational_internal_energy(self, temperature, volume):
        """
        Vibrational Helmholtz free energy, A_vib(V, T).
        Eq(4) in doi.org/10.1016/j.comphy.2003.12.001

        Args:
            temperature (float): temperature in K
            volume (float)

        Returns:
            float: vibrational free energy
        """
        y = self.debye_temperature(volume) / temperature
        return 8.617332 * 1e-5 * self.natoms * temperature * (9. / 8. * y + 3 * self.debye_integral(y))

    def debye_temperature(self, volume):
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
        avg_mass = 1.6605e-27 * self.mass / self.natoms
        term1 = (2. / 3. * (1. + self.poisson) / (1. - 2. * self.poisson)) ** (1.5)
        term2 = (1. / 3. * (1. + self.poisson) / (1. - self.poisson)) ** (1.5)
        f = (3. / (2. * term1 + term2)) ** (1. / 3.)
        return 2.9772e-11 * (volume / self.natoms) ** (-1. / 6.) * f * np.sqrt(self.bulk_modulus/avg_mass)

    def debye_integral(self, y):
        """
        Debye integral. Eq(5) in  doi.org/10.1016/j.comphy.2003.12.001

        Args:
            y (float): debye temperature/T
            integrator: function from scipy.integrate

        Returns:
            float
        """
        # floating point limit is reached around y=155, so values beyond that are set to the
        # limiting value(T-->0, y --> \infty) of 6.4939394 (from wolfram alpha).
        factor = 3. / y ** 3
        if y < 155:
            return list(quadrature(lambda x: x ** 3 / (np.exp(x) - 1.), 0, y))[0] * factor
        else:
            return 6.493939 * factor

    def gruneisen_parameter(self, temperature, volume):
        return volume * self.pressure / self.vibrational_internal_energy(temperature,volume)

    def thermal_conductivity(self, temperature, volume):
        pass

    def summary_dict(self):
        d = defaultdict(list)
        d["gibbs_free_energy"] = self.gibbs_free_energy
        d["temperatures"] = self.temperatures
        d["optimum_volumes"] = self.optimum_volumes
        for v, t in zip(self.optimum_volumes, self.temperatures):
            d["debye_temperature"].append(self.debye_temperature(v))
            d["gruneisen_parameter"].append(self.gruneisen_parameter(t, v))
            d["thermal_conductivity"].append(self.thermal_conductivity(t, v))
        return d


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
