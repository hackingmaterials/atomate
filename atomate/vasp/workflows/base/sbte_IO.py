#!/usr/bin/env python3

import f90nml
"""
This module defines IO ShengBTE for making CONTROL file (i.e. a class that can read/write the ‘control’ file, f90ml)

TO-DO: Test this thing
"""

class ShengBTE_CONTROL_IO:

    """

        Args:
            nelements (int): number of elements in the material
            natoms (int): number of atoms in the cell
            ngrid (3 length list of ints): q-point mesh
            ***norientations (int):
            ***lfactor (float):
            lattvec1 (3 length list of floats): first lattice vector
            lattvec2 (3 length list of floats): second lattice vector
            lattvec3 (3 length list of floats): third lattice vector
            elements (comma-separated list of strings): list of elements present
            types (natoms length list of ints): atom types listed in the same order as in 'elements' (first atom should be '1')
            positions (natoms x 3 matrix): atomic positions of each atom
            scell (3 length list of ints): dimensions of supercell input
            ***temperature (int): temperature in Kelvin?
            scalebroad (float): Gaussian broadening of allowable energies that can scatter a phonon to a certain state
                (higher broadening = more scattering allowed. Too high is unphysical. For complex
                materials, make 'scalebroad' smaller since there will be more scattering candidates
                for any given phonon.
            ***onlyharmonic (bool): Whether to only consider harmonic scattering (2-phonon?)
            ***isotopes (bool): Whether to consider scattering from isotopic disorder.
            ***nonanalytic (bool): Used for Born effective charge and dielectric constant.
            nanowires (bool): Whether the structure input is a nanowire.

        Returns:
            IO reader/writer for ShengBTE's CONTROL file, composed of FORTRAN namelists
    """

    def __init__(self, nelements=None, natoms=None, elements=None, types=None, positions=None, ngrid=[25, 25, 25],
                 lfactor=0.1, lattvec1=[1, 0, 0], lattvec2=[0, 1, 0], lattvec3=[0, 0, 1], scell=[5, 5, 5],
                 temperature=500,  scalebroad=0.5, onlyharmonic=False, isotopes=False, nonanalytic=True,
                 nanowires=False, norientations=0):

        self.nelements = nelements
        self.natoms = natoms
        self.ngrid = ngrid
        self.norientations = norientations
        self.lfactor = lfactor
        self.lattvec1 = lattvec1
        self.lattvec2 = lattvec2
        self.lattvec3 = lattvec3
        self.elements = elements
        self.types = types
        self.positions = positions
        self.scell = scell
        self.temperature = temperature
        self.scalebroad = scalebroad
        self.onlyharmonic = onlyharmonic
        self.isotopes = isotopes
        self.nonanalytic = nonanalytic
        self.nanowires = nanowires


    def read_CONTROL(self, filename):
        nml = f90nml.read(filename)
        return nml

    def write_CONTROL(self):
        with open('CONTROL.nml', 'w') as control_file:
            nml = {'allocations':
                       {'nelements': self.nelements,
                        'natoms': self.natoms,
                        'ngrid': self.ngrid,
                        'norientations': self.norientations},
                   'crystal':
                       {'lfactor': self.lfactor,
                        'lattvec(:,1)': self.lattvec1,
                        'lattvec(:,2)': self.lattvec2,
                        'lattvec(:,3)': self.lattvec3,
                        'elements': self.elements,
                        'types': self.types,
                        'scell': self.scell},
                   'parameters':
                       {'T': self.temperature,
                        'scalebroad': self.scalebroad},
                   'flags':
                       {'onlyharmonic': self.onlyharmonic,
                        'isotopes': self.isotopes,
                        'nonanalytic': self.nonanalytic,
                        'nanowires': self.nanowires}}

            #add positions to the dict
            num_atoms, _ = self.positions.shape
            for at in range(num_atoms):
                at_position = 'positions(:,' + str(at+1) + ')'
                nml['crystal'][at_position] = self.positions[at, :]

            f90nml.write(nml, control_file, force=True) #force=True overwrites an existing file
