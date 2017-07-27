# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines functions that yield lammps workflows
"""

from fireworks import Workflow

from pymatgen.io.lammps.sets import LammpsInputSet

from atomate.lammps.fireworks.core import LammpsFW

__author__ = 'Kiran Mathew'
__email__ = "kmathew@lbl.gov"


def get_wf(name, lammps_input_set, input_filename, data_filename, lammps_cmd, db_file):
    """
    Returns workflow that writes lammps input/data files, runs lammps and inserts to DB.

    Args:
        name (str): workflow name
        lammps_input_set (DictLammpsInput): lammps input set
        input_filename (str): input file name
        data_filename (str): data file name
        lammps_cmd (string): path to the lammps binary
        db_file (string): path to the db file

    Returns:
        Workflow
    """

    fws = [LammpsFW(lammps_input_set=lammps_input_set, input_filename=input_filename,
                    data_filename=data_filename, lammps_cmd=lammps_cmd, db_file=db_file)]
    return Workflow(fws, name=name)


def get_wf_from_input_template(input_template_file, user_settings, lammps_data=None,
                               input_filename="lammps.in", data_filename="lammps.data",
                               is_forcefield=False, lammps_cmd="lammps", db_file=None,
                               name="LAMMPS template Wflow"):
    """
    Returns workflow where the input file parameters are set from the give json(or plain text)
    template file.

    Args:
        input_template_file (str): path to plain text template lammps input file.
        user_settings (dict): User lammps settings
        lammps_data (string/LammpsData/LammpsForceFieldData): path to the data file or
            an appropriate object.
        input_filename (string): input file name. This is the name of the input file passed to the
            lammps binary.
        data_filename (string): name of the the lammps data file.
        is_forcefield (bool): whether the data file has forcefield and topology info in it.
            This is required only if lammps_data is a path to the data file instead of a data object.
        lammps_cmd (string): lammps command to run (skip the input file).
        db_file (string): path to the db file.
        name (str): workflow name

    Returns:
        Workflow
    """
    wf_name = name
    lammps_input_set = LammpsInputSet.from_file(wf_name, input_template_file,
                                                user_settings=user_settings,
                                                lammps_data=lammps_data,
                                                data_filename=data_filename,
                                                is_forcefield=is_forcefield)

    return get_wf(wf_name, lammps_input_set, input_filename=input_filename,
                  data_filename=data_filename, lammps_cmd=lammps_cmd, db_file=db_file)
