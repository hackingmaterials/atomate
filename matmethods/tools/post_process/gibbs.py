"""
Compute the quasi harmonic approximation gibbs free energy using phonopy.
required inputs: volume, energy and the corresponding dynamical matrix.

TODO: cleanup and generalize.
"""

import numpy as np

from phonopy import Phonopy
from phonopy.structure.atoms import Atoms as PhonopyAtoms
from phonopy import PhonopyQHA

import matplotlib.pyplot as plt

__author__ = "Kiran Mathew"
__email__ = "kmathew@lbl.gov"


Si_primitive = PhonopyAtoms(symbols=['Si'] * 2, 
                            scaled_positions=[(0, 0, 0), (0.75, 0.5, 0.75)] )
Si_primitive.set_cell([[3.867422 ,0.000000, 0.000000],
                       [1.933711, 3.349287, 0.000000],
                       [-0.000000, -2.232856, 3.157737]])

# supercell size
scell = [[1,0,0],[0,1,0],[0,0,1]]
phonon = Phonopy(Si_primitive, scell)


eos = "vinet" # birch_murnaghan, and murnaghan
t_step=10
t_min=0
t_max=1000
mesh = [20, 20, 20]

energies =  [-10.74872285, -10.81230037, -10.76735579, -10.80674705, 
             -10.63448195, -10.85052634, -10.54874253]# in eV

volumes  = [ 37.3151026513, 43.3880395077, 44.6767086629, 38.4811269875, 
            47.3300938608, 40.8855198875, 35.0542180638] # in Angs^3

force_constants = [
    
    [[[[-0.58129914, -0.0, -0.0], [-0.0, -0.58129937, -3.8e-07], [-0.0, -3.8e-07, -0.58129926]], [[0.58129914, 0.0, 0.0], [0.0, 0.58129937, 3.8e-07], [0.0, 3.8e-07, 0.58129926]]], [[[0.58129914, 0.0, 0.0], [0.0, 0.58129937, 3.8e-07], [0.0, 3.8e-07, 0.58129926]], [[-0.58129914, -0.0, -0.0], [-0.0, -0.58129937, -3.8e-07], [-0.0, -3.8e-07, -0.58129926]]]],
    
    [[[[-0.43642287, -0.0, -0.0], [-0.0, -0.43642333, -2.7e-07], [-0.0, -2.7e-07, -0.43642303]], [[0.43642287, 0.0, 0.0], [0.0, 0.43642333, 2.7e-07], [0.0, 2.7e-07, 0.43642303]]], [[[0.43642287, 0.0, 0.0], [0.0, 0.43642333, 2.7e-07], [0.0, 2.7e-07, 0.43642303]], [[-0.43642287, -0.0, -0.0], [-0.0, -0.43642333, -2.7e-07], [-0.0, -2.7e-07, -0.43642303]]]],
    
    [[[[-0.41484091, -6e-08, -4e-08], [-6e-08, -0.41484127, -3e-07], [-4e-08, -3e-07, -0.41484103]], [[0.41484091, 6e-08, 4e-08], [6e-08, 0.41484127, 3e-07], [4e-08, 3e-07, 0.41484103]]], [[[0.41484091, 6e-08, 4e-08], [6e-08, 0.41484127, 3e-07], [4e-08, 3e-07, 0.41484103]], [[-0.41484091, -6e-08, -4e-08], [-6e-08, -0.41484127, -3e-07], [-4e-08, -3e-07, -0.41484103]]]],
    
    [[[[-0.55106927, -8e-08, -6e-08], [-8e-08, -0.55106971, -3.6e-07], [-6e-08, -3.6e-07, -0.55106934]], [[0.55106927, 8e-08, 6e-08], [8e-08, 0.55106971, 3.6e-07], [6e-08, 3.6e-07, 0.55106934]]], [[[0.55106927, 8e-08, 6e-08], [8e-08, 0.55106971, 3.6e-07], [6e-08, 3.6e-07, 0.55106934]], [[-0.55106927, -8e-08, -6e-08], [-8e-08, -0.55106971, -3.6e-07], [-6e-08, -3.6e-07, -0.55106934]]]],
    
    [[[[-0.37167229, -5e-08, -4e-08], [-5e-08, -0.37167246, -2.2e-07], [-4e-08, -2.2e-07, -0.37167235]], [[0.37167229, 5e-08, 4e-08], [5e-08, 0.37167246, 2.2e-07], [4e-08, 2.2e-07, 0.37167235]]], [[[0.37167229, 5e-08, 4e-08], [5e-08, 0.37167246, 2.2e-07], [4e-08, 2.2e-07, 0.37167235]], [[-0.37167229, -5e-08, -4e-08], [-5e-08, -0.37167246, -2.2e-07], [-4e-08, -2.2e-07, -0.37167235]]]] ,
    
    [[[[-0.48985069, 7e-08, 5e-08], [7e-08, -0.48985087, -2.6e-07], [5e-08, -2.6e-07, -0.48985069]], [[0.48985069, -7e-08, -5e-08], [-7e-08, 0.48985087, 2.6e-07], [-5e-08, 2.6e-07, 0.48985069]]], [[[0.48985069, -7e-08, -5e-08], [-7e-08, 0.48985087, 2.6e-07], [-5e-08, 2.6e-07, 0.48985069]], [[-0.48985069, 7e-08, 5e-08], [7e-08, -0.48985087, -2.6e-07], [5e-08, -2.6e-07, -0.48985069]]]],
    
    [[[[-0.65665902, 0.0, 0.0], [0.0, -0.65665956, -5.1e-07], [0.0, -5.1e-07, -0.65665925]], [[0.65665902, -0.0, -0.0], [-0.0, 0.65665956, 5.1e-07], [-0.0, 5.1e-07, 0.65665925]]], [[[0.65665902, -0.0, -0.0], [-0.0, 0.65665956, 5.1e-07], [-0.0, 5.1e-07, 0.65665925]], [[-0.65665902, 0.0, 0.0], [0.0, -0.65665956, -5.1e-07], [0.0, -5.1e-07, -0.65665925]]]]
]


# The physical units of V and T are \AA^3 and K, respectively.
# The unit of eV for Helmholtz and Gibbs energies, 
# J/K/mol for C_V and entropy, GPa for for bulk modulus and pressure are used.
temperatures = []
free_energy = []
entropy = []
cv = []

for f in force_constants:
    phonon.set_force_constants(-np.array(force_constants))
    phonon.set_mesh(mesh)
    phonon.set_thermal_properties(t_step=t_step, t_min=t_min, t_max=t_max)
    t, g, e, c = phonon.get_thermal_properties()
    temperatures.append(t)
    free_energy.append(g)
    entropy.append(e)
    cv.append(c)

phonopy_qha = PhonopyQHA(volumes, energies, eos=eos, 
                         temperatures=temperatures[0],
                         free_energy=np.array(free_energy).T, cv=np.array(cv).T,
                         entropy=np.array(entropy).T, 
                         t_max=np.max(temperatures[0]))


# gibbs free energy
max_t_index = phonopy_qha._qha._max_t_index
G = phonopy_qha.get_gibbs_temperature()[:max_t_index]
T = phonopy_qha._qha._temperatures[:max_t_index]

plt.plot(T, G)

phonopy_qha.plot_qha(thin_number=10, volume_temp_exp=None).show()
