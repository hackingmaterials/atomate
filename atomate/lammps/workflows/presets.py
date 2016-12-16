# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines functions that yield lammps workflows
"""

from fireworks import Workflow

from pymatgen.io.lammps.input import DictLammpsInput, NVTLammpsInput

from atomate.lammps.workflows.core import get_wf


__author__ = 'Kiran Mathew'
__email__ = "kmathew@lbl.gov"


def wf_from_input_template(input_template_file, lammps_data, data_filename, user_settings,
                           is_forcefield=False, input_filename="lammps.inp", lammps_bin="lammps",
                           db_file=None, dry_run=False):
    """
    Returns workflow where the input file parameters are set from the give json template file.

    Args:
        input_template_file: json template input file
        lammps_data (string/LammpsData/LammpsForceFieldData): path to the data file or
            an appropriate object
        data_filename (string): name of the the lammps data file
        user_settings (dict): User lammps settings
        is_forcefield (bool): whether the data file has forcefield and topology info in it.
            This is required only if lammps_data is a path to the data file instead of a data object
        input_filename (string): input file name
        lammps_bin (string): path to the lammps binary
        db_file (string): path to the db file
        dry_run (bool): for test purposes, decides whether or not to run the lammps binary
            with the input file.

    Returns:
        Workflow

    """
    wf_name = "LAMMPS Wflow from input template {}".format(input_template_file)
    lammps_dict_input = DictLammpsInput.from_file(wf_name, input_template_file, lammps_data=lammps_data,
                                                  data_filename=data_filename,
                                                  user_lammps_settings=user_settings,
                                                  is_forcefield=is_forcefield)
    return get_wf(wf_name, lammps_dict_input, input_filename=input_filename, lammps_bin=lammps_bin,
                  db_file=db_file, dry_run=dry_run)


def nvt_wf(lammps_data, input_filename = "nvt.inp", data_filename="in.data", user_lammps_settings={},
           is_forcefield=False, lammps_bin="lammps", db_file=None, dry_run=False):
    """
    Returns NVT workflow:
        Firework: [write lammps input task, run direct task]

    Args:
        lammps_data (string/LammpsData/LammpsForceFieldData): path to the data file
            or an appropriate object.
        input_filename (string): input file name
        data_filename (string): data file name
        user_lammps_settings (dict): used to override the default input file
            paramter settings
        is_forcefield (bool): whether or not the data file has forcefiled info.
        lammps_bin (string): path to the lammps binary
        db_file (string): path to the db file
        dry_run (bool): for test purposes, decides whether or not to run the lammps binary
            with the input file.
    """
    wf_name = "LAMMPS NVT"
    lammps_dict_input = NVTLammpsInput(lammps_data=lammps_data, data_filename=data_filename,
                                       user_lammps_settings=user_lammps_settings, is_forcefield=is_forcefield)
    return get_wf(wf_name, lammps_dict_input, input_filename=input_filename, lammps_bin=lammps_bin,
                  db_file=db_file, dry_run=dry_run)
