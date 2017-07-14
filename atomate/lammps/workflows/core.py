# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines functions that yield lammps workflows
"""

from fireworks import Workflow

from pymatgen.io.lammps.input import DictLammpsInput, NVTLammpsInput

from atomate.lammps.fireworks.core import LammpsFW

__author__ = 'Kiran Mathew'
__email__ = "kmathew@lbl.gov"


def get_wf(name, lammps_input_set, input_filename="lammps.in", lammps_cmd="lammps", db_file=None):
    """
    Returns workflow that writes lammps input/data files, runs lammps and inserts to DB.

    Args:
        name (str): workflow name
        lammps_input_set (DictLammpsInput): lammps input set
        input_filename (string): input file name
        lammps_bin (string): path to the lammps binary
        db_file (string): path to the db file

    Returns:
        Workflow

    """

    fws = [LammpsFW(name, lammps_input_set=lammps_input_set, input_filename=input_filename, lammps_cmd=lammps_cmd,
                    db_file=db_file)]
    return Workflow(fws, name=name)


def get_wf_from_input_template(job_name, input_template_file, lammps_data,
                               input_filename="lammps.in", data_filename="lammps.data",
                               user_lammps_settings=None, is_forcefield=False, lammps_cmd="lammps",
                               db_file=None):
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
        lammps_cmd (string): path to the lammps binary
        db_file (string): path to the db file
        dry_run (bool): for test purposes, decides whether or not to run the lammps binary
            with the input file.

    Returns:
        Workflow

    """
    user_lammps_settings = user_lammps_settings or {}
    wf_name = job_name or "LAMMPS Wflow from input template {}".format(input_template_file)
    lammps_dict_input = DictLammpsInput.from_file(wf_name, input_template_file,
                                                  lammps_data=lammps_data,
                                                  data_filename=data_filename,
                                                  user_lammps_settings=user_lammps_settings,
                                                  is_forcefield=is_forcefield)
    return get_wf(wf_name, lammps_dict_input, input_filename=input_filename, lammps_cmd=lammps_cmd,
                  db_file=db_file)


def nvt_wf(job_name, lammps_data, input_filename="nvt.inp", data_filename="in.data",
           user_lammps_settings=None, is_forcefield=False, lammps_cmd="lammps", db_file=None):
    """
    Returns NVT workflow (single Firework: [write lammps input task, run direct task])

    Args:
        lammps_data (string/LammpsData/LammpsForceFieldData): path to the data file
            or an appropriate object.
        input_filename (string): input file name
        data_filename (string): data file name
        user_lammps_settings (dict): used to override the default input file parameter settings
        is_forcefield (bool): whether or not the data file has forcefiled info.
        lammps_bin (string): path to the lammps binary
        db_file (string): path to the db file
        dry_run (bool): for test purposes, decides whether or not to run the lammps binary
            with the input file.
    """
    user_lammps_settings = user_lammps_settings or None
    wf_name = job_name or "LAMMPS NVT"
    lammps_dict_input = NVTLammpsInput(lammps_data=lammps_data, data_filename=data_filename,
                                       user_lammps_settings=user_lammps_settings,
                                       is_forcefield=is_forcefield)
    return get_wf(wf_name, lammps_dict_input, input_filename=input_filename, lammps_cmd=lammps_cmd,
                  db_file=db_file)