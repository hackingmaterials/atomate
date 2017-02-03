# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

from collections import defaultdict

import numpy as np

from scipy.optimize import minimize
from scipy.integrate import quadrature
from scipy.constants import physical_constants
from scipy.misc import derivative

from pymatgen.analysis.eos import EOS

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'


class QuasiharmonicDebyeApprox(object):
    """
    Implements the quasiharmonic Debye approximation as described in papers:
    http://doi.org/10.1016/j.comphy.2003.12.001 (2004) and
    http://doi.org/10.1103/PhysRevB.90.174107 (2014)

    Args:
        energies (list): list of DFT energies in eV
        volumes (list): list of volumes in Ang^3
        structure (Structure):
        t_min (float): min temperature
        t_step (float): temperature step
        t_max (float): max temperature
        eos (str): equation of state used for fitting the energies and the volumes.
            options supported by pymatgen: "quadratic", "murnaghan", "birch", "birch_murnaghan",
            "pourier_tarantola", "vinet", "deltafactor"
        pressure (float): in GPa, optional.
        poisson (float): poisson ratio.
        use_mie_gruneisen (bool): whether or not to use the mie-gruneisen formulation to compute
            the gruneisen parameter. The default is slater-gamma formulation.
    """
    def __init__(self, energies, volumes, structure, t_min=300.0, t_step=100, t_max=300.0,
                 eos="vinet", pressure=0.0, poisson=0.25, use_mie_gruneisen=False):
        self.energies = energies
        self.volumes = volumes
        self.structure = structure
        self.temperature_min = t_min
        self.temperature_max = t_max
        self.temperature_step = t_step
        self.eos = eos
        self.pressure = pressure
        self.poisson = poisson
        self.use_mie_gruneisen = use_mie_gruneisen
        self.mass = sum([e.atomic_mass for e in self.structure.species])
        self.natoms = self.structure.composition.num_atoms
        self.avg_mass = physical_constants["atomic mass constant"][0] * self.mass / self.natoms  # kg
        self.kb = physical_constants["Boltzmann constant in eV/K"][0]
        self.hbar = physical_constants["Planck constant over 2 pi in eV s"][0]
        self.gpa_to_ev_ang = 1./160.21766208  # 1 GPa in ev/Ang^3
        self.gibbs_free_energy = []  # optimized values, eV
        self.temperatures = []  # list of temperatures for which the optimized values are available, K
        self.optimum_volumes = []  # in Ang^3
        # fit E and V and get the bulk modulus(used to compute the debye temperature)
        print("Fitting E and V")
        self.eos = EOS(eos)
        self.ev_eos_fit = self.eos.fit(volumes, energies)
        self.bulk_modulus = self.ev_eos_fit.b0_GPa  # in GPa
        self.optimize_gibbs_free_energy()

    def optimize_gibbs_free_energy(self):
        """
        Evaluate the gibbs free energy as a function of V, T and P i.e G(V, T, P),
        minimize G(V, T, P) wrt V for each T and store the optimum values.

        Note: The data points for which the equation of state fitting fails are skipped.
        """
        temperatures = np.linspace(self.temperature_min,  self.temperature_max,
                                   np.ceil((self.temperature_max-self.temperature_min)/self.temperature_step)+1)
        print("Fitting G and V for each T")
        for t in temperatures:
            try:
                G_opt, V_opt = self.optimizer(t)
            except:
                print("EOS fitting failed, so skipping this data point")
                continue
            self.gibbs_free_energy.append(G_opt)
            self.temperatures.append(t)
            self.optimum_volumes.append(V_opt)

    def optimizer(self, temperature):
        """
        Evaluate G(V, T, P) at the given temperature(and pressure) and minimize it wrt V.

        1. Compute the  vibrational helmholtz free energy, A_vib.
        2. Compute the gibbs free energy as a function of volume, temperature and pressure, G(V,T,P).
        3. Preform an equation of state fit to get the functional form of gibbs free energy:G(V, T, P).
        4. Finally G(V, P, T) is minimized with respect to V.

        Args:
            temperature (float): temperature in K

        Returns:
            float, float: G_opt(V_opt, T, P) in eV and V_opt in Ang^3.
        """
        G_V = []  # G for each volume
        # G = E(V) + PV + A_vib(V, T)
        for i, v in enumerate(self.volumes):
            G_V.append(self.energies[i] +
                       self.pressure * v * self.gpa_to_ev_ang +
                       self.vibrational_free_energy(temperature, v))

        # fit equation of state, G(V, T, P)
        eos_fit = self.eos.fit(self.volumes, G_V)
        params = tuple(eos_fit.eos_params.tolist())  # E0(ref energy), B0, B1, V0(ref volume)
        # minimize the fit eos wrt volume
        # Note: the ref energy and the ref volume(E0 and V0) not necessarily the same as
        # minimum energy and min volume.
        min_wrt_vol = minimize(eos_fit.func, min(eos_fit.volumes), params)
        # G_opt=G(V_opt, T, P), V_opt
        return min_wrt_vol.fun, min_wrt_vol.x[0]

    def vibrational_free_energy(self, temperature, volume):
        """
        Vibrational Helmholtz free energy, A_vib(V, T).
        Eq(4) in doi.org/10.1016/j.comphy.2003.12.001

        Args:
            temperature (float): temperature in K
            volume (float)

        Returns:
            float: vibrational free energy in eV
        """
        y = self.debye_temperature(volume) / temperature
        return self.kb * self.natoms * temperature * (9. / 8. * y + 3 * np.log(1 - np.exp(-y))
                                                      - self.debye_integral(y))

    def vibrational_internal_energy(self, temperature, volume):
        """
        Vibrational internal energy, U_vib(V, T).
        Eq(4) in doi.org/10.1016/j.comphy.2003.12.001

        Args:
            temperature (float): temperature in K
            volume (float): in Ang^3

        Returns:
            float: vibrational internal energy in eV
        """
        y = self.debye_temperature(volume) / temperature
        return self.kb * self.natoms * temperature * (9. / 8. * y + 3 * self.debye_integral(y))

    def debye_temperature(self, volume):
        """
        Calculates the debye temperature.
        Eq(6) in doi.org/10.1016/j.comphy.2003.12.001. Thanks to Joey.

        Args:
            volume (float): in Ang^3

        Returns:
            float: debye temperature in K
         """
        term1 = (2. / 3. * (1. + self.poisson) / (1. - 2. * self.poisson)) ** (1.5)
        term2 = (1. / 3. * (1. + self.poisson) / (1. - self.poisson)) ** (1.5)
        f = (3. / (2. * term1 + term2)) ** (1. / 3.)
        return 2.9772e-11 * (volume / self.natoms) ** (-1. / 6.) * f * np.sqrt(self.bulk_modulus/self.avg_mass)

    @staticmethod
    def debye_integral(y):
        """
        Debye integral. Eq(5) in  doi.org/10.1016/j.comphy.2003.12.001

        Args:
            y (float): debye temperature/T, upper limit

        Returns:
            float: unitless
        """
        # floating point limit is reached around y=155, so values beyond that are set to the
        # limiting value(T-->0, y --> \infty) of 6.4939394 (from wolfram alpha).
        factor = 3. / y ** 3
        if y < 155:
            return list(quadrature(lambda x: x ** 3 / (np.exp(x) - 1.), 0, y))[0] * factor
        else:
            return 6.493939 * factor

    def gruneisen_parameter(self, temperature, volume):
        """
        Slater-gamma formulation(the default):
            gruneisen paramter = - d log(theta)/ d log(V)
                               = - ( 1/6 + 0.5 d log(B)/ d log(V) )
                               = - (1/6 + 0.5 V/B dB/dV), where dB/dV = d^2E/dV^2 + V * d^3E/dV^3

        Mie-gruneisen formulation:
            Eq(31) in doi.org/10.1016/j.comphy.2003.12.001
            Eq(7) in MA Blanc0 ef al.lJoumal of Molecular Structure (Theochem) 368 (1996) 245-255
            Also se J.-P. Poirier, Introduction to the Physics of the Earth’s Interior, 2nd ed.
                (Cambridge University Press, Cambridge, 2000) Eq(3.53)

        Args:
            temperature (float): temperature in K
            volume (float): in Ang^3

        Returns:
            float: unitless
        """
        # Mie-gruneisen formulation
        if self.use_mie_gruneisen:
            # first derivative of energy at 0K wrt volume evaluated at the given volume, in eV/Ang^3
            p0 = derivative(self.ev_eos_fit.func, volume, dx=1e-3,
                            args=tuple(tuple(self.ev_eos_fit.eos_params.tolist())))
            return (self.gpa_to_ev_ang * volume * (self.pressure + p0 / self.gpa_to_ev_ang) /
                    self.vibrational_internal_energy(temperature, volume))

        # Slater-gamma formulation
        # second derivative of energy at 0K wrt volume evaluated at the given volume, in eV/Ang^6
        d2EdV2 = derivative(self.ev_eos_fit.func, volume, dx=1e-3, n=2,
                            args=tuple(tuple(self.ev_eos_fit.eos_params.tolist())), order=5)
        # third derivative of energy at 0K wrt volume evaluated at the given volume, in eV/Ang^9
        d3EdV3 = derivative(self.ev_eos_fit.func, volume, dx=1e-3, n=3,
                            args=tuple(tuple(self.ev_eos_fit.eos_params.tolist())), order=7)
        # first derivative of bulk modulus wrt volume, eV/Ang^6
        dBdV = d2EdV2 + d3EdV3 * volume
        return -(1./6. + 0.5 * volume * dBdV / self.ev_eos_fit.b0)

    def thermal_conductivity(self, temperature, volume):
        """
        Eq(17) in 10.1103/PhysRevB.90.174107

        Args:
            temperature (float): temperature in K
            volume (float): in Ang^3

        Returns:
            float: thermal conductivity in W/K/m
        """
        gamma = self.gruneisen_parameter(temperature, volume)
        theta_d = self.debye_temperature(volume)  # K
        theta_a = theta_d * self.natoms**(-1./3.)  # K
        prefactor = (0.849 * 3 * 4**(1./3.)) / (20. * np.pi**3)
        prefactor = prefactor * (self.kb/self.hbar)**3 * self.avg_mass  # kg/K^3/s^3
        kappa = prefactor / (gamma**2 - 0.514 * gamma + 0.228)
        # kg/K/s^3 * Ang = (kg m/s^2)/(Ks)*1e-10 = N/(Ks)*1e-10 = Nm/(Kms)*1e-10 = W/K/m*1e-10
        kappa = kappa * theta_a**2 * volume**(1./3.) * 1e-10
        return kappa

    def get_summary_dict(self):
        """
        Returns a dict with a summary of the computed properties.
        """
        d = defaultdict(list)
        d["pressure"] = self.pressure
        d["poisson"] = self.poisson
        d["mass"] = self.mass
        d["natoms"] = int(self.natoms)
        d["bulk_modulus"] = self.bulk_modulus
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
