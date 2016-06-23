# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

from fireworks import Workflow, Firework

from pymatgen.io.lammps.input import NVTLammpsInput

from matmethods.lammps.firetasks.write_inputs import WritelammpsInputFromDictInput
from matmethods.lammps.firetasks.run_calc import RunLammpsDirect


__author__ = 'Kiran Mathew'


def nvt_wf(data_file = "nvt.data", input_file = "nvt.inp"):
    lammps_dict_input = NVTLammpsInput(data_file=data_file, is_forcefield=True)
    task1 = WritelammpsInputFromDictInput(lammps_dict_input=lammps_dict_input, input_file=input_file)
    task2 = RunLammpsDirect(lammps_cmd="lammps -in "+input_file)
    fw1 = Firework([task1, task2], name='Run lammps')
    return Workflow([fw1], name="lammps nvt")
