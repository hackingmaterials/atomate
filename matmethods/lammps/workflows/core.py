# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

from fireworks import Workflow, Firework

from pymatgen.io.lammps.input import NVTLammpsInput

from matmethods.lammps.firetasks.write_inputs import WritelammpsInputFromDictInput
from matmethods.lammps.firetasks.run_calc import RunLammpsDirect


__author__ = 'Kiran Mathew'


def nvt_wf(data_input, input_filename = "nvt.inp", data_filename="in.data",
           user_lammps_settings={}, is_forcefield=False, lammps_bin="lammps"):
    """
    Returns NVT workflow

    Args:
        data_input (string/LammpsData/LammpsForceFieldData): path to the data file
            or an appropriate object.
        input_filename (string): input file name
        data_filename (string): data file name
        user_lammps_settings (dict): used to override the default input file
            paramter settings
        is_forcefield (bool): whether or not the data file has forcefiled info.
        lammps_bin (string): path to the lammps binary
    """
    lammps_dict_input = NVTLammpsInput(lammps_data=data_input, data_filename=data_filename,
                                       user_lammps_settings=user_lammps_settings, is_forcefield=is_forcefield)
    task1 = WritelammpsInputFromDictInput(lammps_dict_input=lammps_dict_input, input_file=input_filename)
    task2 = RunLammpsDirect(lammps_cmd=lammps_bin+" -in "+input_filename)
    fw1 = Firework([task1, task2], name='Run lammps')
    return Workflow([fw1], name="LAMMPS NVT")


if __name__ == "__main__":
    wf = nvt_wf("test_files/nvt.data", data_filename="nvt.data", is_forcefield=True, lammps_bin="lmp_serial")
    print(wf.as_dict())
