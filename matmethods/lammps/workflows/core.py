# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines functions that yield lammps workflows
"""

from fireworks import Workflow, Firework

from pymatgen.io.lammps.input import DictLammpsInput

from matmethods.lammps.firetasks.write_inputs import WritelammpsInputFromDictInput
from matmethods.lammps.firetasks.run_calc import RunLammpsDirect
from matmethods.lammps.firetasks.parse_outputs import  LammpsToDBTask


__author__ = 'Kiran Mathew'
__email__ = "kmathew@lbl.gov"


def get_wf(job_name, lammps_dict_input, input_filename="lammps.inp", lammps_bin="lammps"):
    """
    Returns workflow that writes lammps input/data files, runs lammps and inserts to DB.

    Args:
        job_name: job name
        lammps_dict_input (DictLammpsInput): lammps input
        input_filename (string): input file name
        lammps_bin (string): path to the lammps binary

    Returns:
        Workflow

    """
    task1 = WritelammpsInputFromDictInput(lammps_dict_input=lammps_dict_input, input_file=input_filename)
    task2 = RunLammpsDirect(lammps_cmd=lammps_bin + " -in " + input_filename)
    task3 = LammpsToDBTask(lammps_input=lammps_dict_input)
    fw1 = Firework([task1, task2, task3], name=job_name)
    return Workflow([fw1], name=job_name)
