# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import datetime

from pymatgen.core.structure import Structure
from atomate.utils.utils import get_logger
from atomate import __version__ as atomate_version

logger = get_logger(__name__)

class McsqsDrone:

    def assimilate(self, path):

        logger.info("Assimilating mcsqs for base dir: {}".format(path))

        input_structure = Structure.from_file('rndstr.in')
        output_structure = Structure.from_file('bestsqs.out')

        # get our objective function
        with open('bestcorr.out') as f:
            lines = f.read().split('\n')
            lines = [l for l in lines if l]
            try:
                obj_func = float(lines[-1].split('=')[1])
            except:
                obj_func = str(lines[-1].split('=')[1])

        # get scaling matrix
        try:
            l = input_structure.lattice
            l_sqs = output_structure.lattice
            scaling_matrix = l.find_mapping(l_sqs)[2].tolist()
        except:
            scaling_matrix = "Could not determine scaling matrix?"

        # get mcsqs version
        try:
            # if running via atomate, we will create this file
            with open('mcsqs_version.txt') as f:
                mcsqs_version = str(f.read())
        except:
            mcsqs_version = "Unknown mcsqs version"

        # get number of clusters
        with open('clusters.out') as f:
            lines = f.read().split('\n')
            num_clusters = (len(lines.split('\n')) - 2) / 7

        return {
            'disordered': input_structure,
            'bestsqs': output_structure,
            'clusters': None, # TODO: retrieve this from FW.json
            'num_clusters': num_clusters,
            'objective_function': obj_func,
            'walltime': None, # TODO: pick this up from FW
            'atomate_version': atomate_version,
            'mcsqs_version': mcsqs_version,
            'spacegroup': input_structure.get_space_group_info()[0],
            'scaling_matrix': scaling_matrix,
            'size': len(output_structure)/len(input_structure),
            'last_updated': str(datetime.datetime.now())
        }