# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

from fireworks import FireTaskBase, explicit_serialize


__author__ = 'Kiran Mathew'


@explicit_serialize
class WritelammpsInputFromDictInput(FireTaskBase):
    """
    Writes LAMMPS Input files from DictLammpsInput

    required_params:
        lammps_dict_input (DictLammpsInput)
        input_file (string): path to the input file

    optional_params:
        data_file (string): if specified the data file will be renamed
    """

    required_params = ["lammps_dict_input", "input_file"]
    optional_params = ["data_file"]

    def run_task(self, fw_spec):
        lammps_input = self["lammps_dict_input"]
        lammps_input.write_input(self["input_file"], data_file=self.get("data_file", None))


