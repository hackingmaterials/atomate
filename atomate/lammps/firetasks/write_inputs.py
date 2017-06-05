# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines firetasks for writing LAMMPS input files (data file and the control 
parameters file)
"""

from fireworks import FiretaskBase, explicit_serialize


__author__ = 'Kiran Mathew'
__email__ = "kmathew@lbl.gov"


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

    required_params = ["lammps_input_set", "input_file"]
    optional_params = ["data_file"]

    def run_task(self, fw_spec):
        lammps_input = self["lammps_input_set"]
        lammps_input.write_input(self["input_file"], data_filename=self.get("data_file", None))
