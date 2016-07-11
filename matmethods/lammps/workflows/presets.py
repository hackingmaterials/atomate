# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines functions that yield lammps workflows
"""

from fireworks import Workflow

from pymatgen.io.lammps.input import DictLammpsInput, NVTLammpsInput

from matmethods.lammps.workflows.core import get_wf


__author__ = 'Kiran Mathew'
__email__ = "kmathew@lbl.gov"


def wf_from_input_template(input_template_file, lammps_data, data_filename, user_settings,
                           is_forcefield=False, input_filename="lammps.inp", lammps_bin="lammps"):
    """
    Returns workflow where the input file paramters are set from the give json template file.

    Args:
        job_name: job name
        input_template_file: json template file
        lammps_data (string/LammpsData/LammpsForceFieldData): path to the
                data file or an appropriate object
        data_filename (string): name of the the lammps data file
        user_settings (dict): User lammps settings
        is_forcefield (bool): whether the data file has forcefield and
                topology info in it. This is required only if lammps_data is
                a path to the data file instead of a data object
        input_filename (string): input file name
        lammps_bin (string): path to the lammps binary

    Returns:
        Workflow

    """
    wf_name = "LAMMPS Wflow from input template {}".format(input_template_file)
    lammps_dict_input = DictLammpsInput.from_file(job_name, input_template_file, lammps_data=lammps_data,
                                                  data_filename=data_filename,
                                                  user_lammps_settings=user_settings,
                                                  is_forcefield=is_forcefield)
    return get_wf(wf_name, lammps_dict_input, input_filename=input_filename, lammps_bin=lammps_bin)


def nvt_wf(data_input, input_filename = "nvt.inp", data_filename="in.data",
           user_lammps_settings={}, is_forcefield=False, lammps_bin="lammps"):
    """
    Returns NVT workflow:
        Firework: [write lammps input task, run direct task]

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
    wf_name = "LAMMPS NVT"
    lammps_dict_input = NVTLammpsInput(lammps_data=data_input, data_filename=data_filename,
                                       user_lammps_settings=user_lammps_settings, is_forcefield=is_forcefield)
    return get_wf(wf_name, lammps_dict_input, input_filename=input_filename, lammps_bin=lammps_bin)


if __name__ == "__main__":
    wf = nvt_wf("test_files/nvt.data", data_filename="nvt.data", is_forcefield=True, lammps_bin="lmp_serial")
    print(wf.as_dict())
