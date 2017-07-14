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

    required_params = ["job_name", "lammps_input",  "lammps_data", "is_forcefield"]

    optional_params = ["user_lammps_settings", "data_filename"]

    def run_task(self, fw_spec):

        lammps_input = self["lammps_input"]
        lammps_data = self["lammps_data"]
        job_name = self["job_name"]
        is_forcefield = self["is_forcefield"]
        user_lammps_settings = self.get("user_lammps_settings", {"log": "lammps.log"})
        data_filename = self.get("data_filename", "lammps.data")

        lammps_input_set = DictLammpsInput.from_file(job_name, lammps_input, lammps_data=lammps_data,
                                                     data_filename=data_filename,
                                                     user_lammps_settings=user_lammps_settings,
                                                     is_forcefield=is_forcefield)
        lammps_input_set.write_input("lammps.in", data_filename)
