# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines tasks for writing vasp input sets for various types of vasp calculations
"""

import os
from six.moves import range
from importlib import import_module

import numpy as np

from monty.serialization import dumpfn

from fireworks import FiretaskBase, explicit_serialize

from atomate.utils.utils import env_chk, load_class

__author__ = 'Nicholas Winner'
__email__ = 'ajain@lbl.gov'


@explicit_serialize
class WriteCp2kFromIOSet(FiretaskBase):
    """
    Create VASP input files using implementations of pymatgen's AbstractVaspInputSet. An input set
    can be provided as an object or as a String/parameter combo.

    Required params:
        structure (Structure): structure
        vasp_input_set (AbstractVaspInputSet or str): Either a VaspInputSet object or a string
            name for the VASP input set (e.g., "MPRelaxSet").

    Optional params:
        vasp_input_params (dict): When using a string name for VASP input set, use this as a dict
            to specify kwargs for instantiating the input set parameters. For example, if you want
            to change the user_incar_settings, you should provide: {"user_incar_settings": ...}.
            This setting is ignored if you provide the full object representation of a VaspInputSet
            rather than a String.
    """

    required_params = ["structure", "cp2k_input_set"]
    optional_params = ["cp2k_input_params"]

    def run_task(self, fw_spec):

        if isinstance(self['cp2k_input_set'], dict):
            cis = load_class('pymatgen.io.cp2k.sets',
                             self['cp2k_input_set']['@module']).from_dict(
                             self['cp2k_input_set'])

        # if VaspInputSet String + parameters was provided
        else:
            cis_cls = load_class("pymatgen.io.cp2k.sets", self["cp2k_input_set"])
            cis = cis_cls(self["structure"], **self.get("cp2k_input_params", {}))
        cis.write_file(input_filename='cp2k.inp', output_dir='.')

