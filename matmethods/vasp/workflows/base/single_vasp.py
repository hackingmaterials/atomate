# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines functions that generate workflows for single structural optimization
calculation.
"""

from pymatgen.io.vasp.sets import MPVaspInputSet

from fireworks import Firework, Workflow

from matmethods.common.firetasks.parse_outputs import ToDbTask
from matmethods.vasp.firetasks.run_calc import RunVaspDirect
from matmethods.vasp.firetasks.write_inputs import WriteVaspFromIOSet
from matmethods.vasp.drones import VaspDrone

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


def get_wf_single(structure, vasp_input_set=None, vasp_cmd="vasp", db_file=None,
                  task_label="single VASP"):
    """
    Return vasp workflow consisting of a single firework made of 3 firetasks:
        write inputset for structural relaxation,
        run vasp
        database insertion

    Args:
        structure (Structure): input structure to be relaxed.
        vasp_input_set (DictVaspInputSet): vasp input set.
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.
        task_label (str): workflow name.

    Returns:
        Workflow
    """
    vasp_input_set = vasp_input_set or MPVaspInputSet(force_gamma=True)

    write_task = WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set)
    run_task = RunVaspDirect(vasp_cmd=vasp_cmd)
    parse_task = ToDbTask(db_file=db_file, drone=VaspDrone())

    my_fw = Firework([write_task, run_task, parse_task],
                     name="{}-{}".format(structure.composition.reduced_formula, task_label))
    my_wf = Workflow.from_Firework(my_fw)

    return my_wf

if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest
    structure = PymatgenTest.get_structure("Si")
    wf = get_wf_single(structure)