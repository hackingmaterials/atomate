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

    def read_CONTROL(self, filename):
        nml = f90nml.read(filename)
        return nml

    def write_CONTROL_from_dict(self, dict, filename='CONTROL'):
        with open(filename, 'w') as control_file:
            nml = {'allocations':
                       {'nelements': dict.get('allocations', None).get('nelements', None),
                        'natoms': dict.get('allocations', None).get('natoms', None),
                        'ngrid': dict.get('allocations', [25, 25, 25]).get('ngrid', [25, 25, 25]),
                        'norientations': dict.get('allocations', 0).get('norientations', 0)},
                   'crystal':
                       {'lfactor': dict.get('crystal', 0.1).get('lfactor', 0.1),
                        # 'lattvec(:,1)': dict.get('crystal', None).get('lattvec(:,1)', None),
                        # 'lattvec(:,2)': dict.get('crystal', None).get('lattvec(:,2)', None),
                        # 'lattvec(:,3)': dict.get('crystal', None).get('lattvec(:,3)', None),
                        'lattvec': dict.get('crystal', None).get('lattvec', None),
                        'positions': dict.get('crystal', None).get('positions', None),
                        'elements': dict.get('crystal', None).get('elements', None),
                        'types': dict.get('crystal', None).get('types', None),
                        'scell': dict.get('crystal', [5, 5, 5]).get('scell', [5, 5, 5])},
                   'parameters':
                       {'T': dict.get('parameters', 500).get('T', 500),
                        'scalebroad': dict.get('parameters', 0.5).get('scalebroad', 0.5)},
                   'flags':
                       {'onlyharmonic': dict.get('flags', False).get('onlyharmonic', False),
                        'isotopes': dict.get('flags', False).get('isotopes', False),
                        'nonanalytic': dict.get('flags', True).get('nonanalytic', True),
                        'nanowires': dict.get('flags', False).get('nanowires', False)}}

            f90nml.write(nml, control_file, force=True) #force=True overwrites an existing file


def main():
    # Read the CONTROL file into a Fortran namelist object
    sbte_io = ShengBTE_CONTROL_IO()
    namelist = sbte_io.read_CONTROL('CONTROL')
    print(namelist)
    print('---------------------')

    # Convert the namelist object to a dict for easy access of contents
    dict = namelist.todict()
    print(dict)
    print(dict['allocations']['nelements'])

    # Write the dict back into a namelist file
    sbte_io.write_CONTROL_from_dict(dict, filename='CONTROL_test')


if __name__ == '__main__':
    main()