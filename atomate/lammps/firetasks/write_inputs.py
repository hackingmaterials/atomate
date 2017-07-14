# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines firetasks for writing LAMMPS input files (data file and the control 
parameters file)
"""

from fireworks import FiretaskBase, explicit_serialize
from pymatgen.io.lammps.input import DictLammpsInput


__author__ = 'Kiran Mathew, Brandon Wood'
__email__ = "kmathew@lbl.gov, b.wood@berkeley.edu"


@explicit_serialize
class WriteLammpsFromIOSet(FiretaskBase):
    """
    Writes LAMMPS Input files(data file and the control parameters file) from DictLammpsInput.

    required_params:
        lammps_input_set (DictLammpsInput)
        input_file (string): name of the file to which the input params will be written

    optional_params:
        data_file (string): if specified the data file will be renamed
    """

    required_params = ["lammps_input_set"]

    optional_params = ["input_filename"]

    def run_task(self, fw_spec):

        lammps_input_set = self["lammps_input_set"]
        input_filename = self.get("input_filename", "lammps.in")
        data_filename = self.get("data_filename", "lammps.data")

        lammps_input_set.write_input(input_filename, data_filename)
