# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines tasks for writing vasp input sets for various types of vasp calculations
"""

import os
from six.moves import range
from importlib import import_module
import glob

import numpy as np

from monty.serialization import dumpfn

from pymatgen.io.cp2k.sets import Cp2kInputSet
from pymatgen.io.cp2k.outputs import Cp2kOutput
from fireworks import FiretaskBase, explicit_serialize

from atomate.utils.utils import env_chk, load_class
from atomate.common.firetasks.glue_tasks import get_calc_loc

__author__ = 'Nicholas Winner'
__email__ = 'nwinner@berkeley.edu'


@explicit_serialize
class WriteCp2kFromIOSet(FiretaskBase):
    """
    Create CP2K input files using the pymatgen Cp2kInputSet object or dict representation.

    Required params:
        structure (Structure): structure
        cp2k_input_set (Cp2kInputSet or dict): Either a Cp2kInputSet object or a dict representation of one

    Optional params:
        cp2k_input_params (dict): When using a string name for CP2K input set, use this as a dict
            to specify kwargs for instantiating the input set parameters. For example, if you wan
    """

    required_params = ["structure", "cp2k_input_set"]
    optional_params = ["cp2k_input_params"]

    def run_task(self, fw_spec):
        if isinstance(self['cp2k_input_set'], dict):
            cis = load_class('pymatgen.io.cp2k.sets',
                             self['cp2k_input_set']['@module']).from_dict(
                             self['cp2k_input_set'])
        else:
            cis_cls = load_class("pymatgen.io.cp2k.sets", self["cp2k_input_set"])
            cis = cis_cls(self["structure"], **self.get("cp2k_input_params", {}))
        cis.write_file(input_filename='cp2k.inp', output_dir='.')


@explicit_serialize
class WriteCp2kFromPrevious(FiretaskBase):

    optional_params = ['cp2k_input_params', 'prev_calc_loc', 'original_input_filename',
                       'new_input_filename']

    def run_task(self, fw_spec):
        calc_loc = get_calc_loc('prev_calc_loc', fw_spec["calc_locs"])
        input_filename = self.get('original_input_filename', 'cp2k.input')
        if os.path.isfile(calc_loc, input_filename):
            input_path = os.path.join(calc_loc, input_filename)
        else:
            raise FileNotFoundError("Could not find the cp2k input file!")

        ci = Cp2kInputSet.from_file(input_path)
        ci.write_file(input_filename='cp2k.inp')




