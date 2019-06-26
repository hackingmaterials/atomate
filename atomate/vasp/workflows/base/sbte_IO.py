#!/usr/bin/env python3

import f90nml
import numpy as np
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



    def file_writer_helper_func(self, dict, filename):
        nelements = str(dict['allocations']['nelements'])
        natoms = str(dict['allocations']['natoms'])
        ngrid = dict['allocations']['ngrid']
        norientations = str(dict['allocations']['norientations'])

        lfactor = str(dict['crystal']['lfactor'])
        lattvec1 = dict['crystal']['lattvec'][0]
        lattvec2 = dict['crystal']['lattvec'][1]
        lattvec3 = dict['crystal']['lattvec'][2]
        elements = dict['crystal']['elements'] #string or list of strings, needs to be parsed
        types = dict['crystal']['types'] #list of ints, needs to be parsed
        positions = np.asarray(dict['crystal']['positions']) #numpy array, needs to be parsed
        num_sites, _ = positions.shape
        scell = dict['crystal']['scell'] #list of ints, needs to be parsed

        temperature = str(int(dict['parameters']['T']))
        scalebroad = str(dict['parameters']['scalebroad'])

        onlyharmonic = dict['flags']['onlyharmonic']
        isotopes = dict['flags']['isotopes']
        nonanalytic = dict['flags']['nonanalytic']
        nanowires = dict['flags']['nanowires']

        def boolean_to_string(boolean):
            if boolean is True:
                return '.TRUE.'
            else:
                return '.FALSE.'

        indent = '        '
        types_string ='types='
        positions_string = ''
        for line in range(num_sites):
            if line != num_sites-1:
                types_string += str(types[line])+' '
            else:
                types_string += str(types[line])+',\n'
            positions_string += 'positions(:,' + str(line+1) + ')=' + str(positions[line,0]) + '  ' \
                                + str(positions[line,1]) + '  ' + str(positions[line,2]) + ',\n' + indent

        full_string = '&allocations\n'+indent+'nelements='+nelements+',\n'
        full_string += indent+'natoms='+natoms+',\n'
        full_string += indent+'ngrid(:)='+str(ngrid[0])+' '+str(ngrid[1])+' '+str(ngrid[2])+'\n'
        full_string += indent+'norientations='+norientations+'\n'
        full_string += '&end\n&crystal\n'
        full_string += indent+'lfactor='+lfactor+',\n'
        full_string += indent+'lattvec(:,1)='+str(lattvec1[0])+'  '+str(lattvec1[1])+'  '+str(lattvec1[2])+',\n'
        full_string += indent+'lattvec(:,2)='+str(lattvec2[0])+'  '+str(lattvec2[1])+'  '+str(lattvec2[2])+',\n'
        full_string += indent+'lattvec(:,3)='+str(lattvec3[0])+'  '+str(lattvec3[1])+'  '+str(lattvec3[2])+',\n'
        full_string += indent+'elements='
        if isinstance(elements, list):
            for i in range(len(elements)):
                full_string += '\"'+elements[i]+str('\"')
                if i != (len(elements)-1):
                    full_string += ' '
                else:
                    full_string += '\n'
        else:
            full_string += '\"'+elements+str('\"\n')
        full_string += indent+types_string
        full_string += indent+positions_string
        full_string += 'scell(:)='+str(scell[0])+' '+str(scell[1])+' '+str(scell[2])+'\n'
        full_string += '&end\n&parameters\n'
        full_string += indent+'T='+temperature+'\n'
        full_string += indent+'scalebroad='+scalebroad+'\n'
        full_string += '&end\n&flags\n'
        full_string += indent+'isotopes='+boolean_to_string(isotopes)+'\n'
        full_string += indent+'onlyharmonic='+boolean_to_string(onlyharmonic)+'\n'
        full_string += indent+'nonanalytic='+boolean_to_string(nonanalytic)+'\n'
        full_string += indent+'nanowires='+boolean_to_string(nanowires)+'\n'
        full_string += '&end'

        file = open(filename, 'w+')
        file.write(full_string)
        file.close()


    def write_CONTROL_from_dict(self, dict, filename='CONTROL', overwrite=True):
        with open(filename, 'w') as control_file:
            new_dict = {'allocations':
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
                                'elements': dict.get('crystal', None).get('elements', None),
                                'types': dict.get('crystal', None).get('types', None),
                                'positions': dict.get('crystal', None).get('positions', None),
                                'scell': dict.get('crystal', [5, 5, 5]).get('scell', [5, 5, 5])},
                        'parameters':
                               {'T': dict.get('parameters', 300).get('t', 300),
                                'scalebroad': dict.get('parameters', 0.5).get('scalebroad', 0.5)},
                        'flags':
                               {'onlyharmonic': dict.get('flags', False).get('onlyharmonic', False),
                                'isotopes': dict.get('flags', False).get('isotopes', False),
                                'nonanalytic': dict.get('flags', True).get('nonanalytic', True),
                                'nanowires': dict.get('flags', False).get('nanowires', False)}}

            # nml = f90nml.namelist.Namelist(new_dict) #convert dict to namelist object
            # nml.false_repr = '.FALSE.'
            # nml.true_repr = '.TRUE.'
            # nml.index_spacing = False
            # nml.uppercase = False
            # nml.indent = '        '
            # f90nml.write(nml, control_file, force=overwrite, sort=False) #force=True overwrites an existing file
            self.file_writer_helper_func(new_dict, filename=filename)


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
    print(dict['crystal']['lattvec'])
    print(dict['crystal']['lattvec'][0])

    # Write the dict back into a namelist file
    sbte_io.write_CONTROL_from_dict(dict, filename='CONTROL_test2')


if __name__ == '__main__':
    main()