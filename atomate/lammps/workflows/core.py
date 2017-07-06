# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines functions that yield lammps workflows
"""

from fireworks import Workflow, Firework

from pymatgen.io.lammps.input import DictLammpsInput, NVTLammpsInput

from atomate.lammps.firetasks.write_inputs import WriteLammpsFromIOSet
from atomate.lammps.firetasks.run_calc import RunLammpsDirect
from atomate.lammps.firetasks.parse_outputs import  LammpsToDBTask


__author__ = 'Kiran Mathew'
__email__ = "kmathew@lbl.gov"


# TODO: @matk86 - It is better if you can have a Fireworks module that includes packmol and
# lammps Fireworks. I understand the firework definitions are simple, but it will be easier for
# people to find the Fireworks there since it will be a familiar subpackage structure.
# It should be an easy mod and shouldn't get in the way much. -computron
# TODO: @matk86 - is there any workflow or firework taking into account Packmol?  -computron

def get_wf(job_name, lammps_input_set, input_filename="lammps.inp", lammps_bin="lammps",
           db_file=None, dry_run=False):
    """
    Returns workflow that writes lammps input/data files, runs lammps and inserts to DB.

    Args:
        job_name: job name
        lammps_input_set (DictLammpsInput): lammps input
        input_filename (string): input file name
        lammps_bin (string): path to the lammps binary
        db_file (string): path to the db file
        dry_run (bool): for test purposes, decides whether or not to run the lammps binary
            with the input file.

    Returns:
        Workflow

    """
    task1 = WriteLammpsFromIOSet(lammps_input_set=lammps_input_set, input_file=input_filename)
    if dry_run:
        lammps_cmd = lammps_bin
    else:
        lammps_cmd = lammps_bin + " -in " + input_filename
    task2 = RunLammpsDirect(lammps_cmd=lammps_cmd)
    task3 = LammpsToDBTask(lammps_input=lammps_input_set, db_file=db_file)
    fw1 = Firework([task1, task2, task3], name=job_name)
    return Workflow([fw1], name=job_name)


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
    lammps_dict_input = DictLammpsInput.from_file(wf_name, input_template_file,
                                                  lammps_data=lammps_data,
                                                  data_filename=data_filename,
                                                  user_lammps_settings=user_settings,
                                                  is_forcefield=is_forcefield)
    return get_wf(wf_name, lammps_dict_input, input_filename=input_filename, lammps_bin=lammps_bin,
                  db_file=db_file, dry_run=dry_run)


def nvt_wf(lammps_data, input_filename = "nvt.inp", data_filename="in.data",
           user_lammps_settings={}, is_forcefield=False, lammps_bin="lammps", db_file=None,
           dry_run=False):
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
    wf_name = "LAMMPS NVT"
    lammps_dict_input = NVTLammpsInput(lammps_data=lammps_data, data_filename=data_filename,
                                       user_lammps_settings=user_lammps_settings,
                                       is_forcefield=is_forcefield)
    return get_wf(wf_name, lammps_dict_input, input_filename=input_filename, lammps_bin=lammps_bin,
                  db_file=db_file, dry_run=dry_run)