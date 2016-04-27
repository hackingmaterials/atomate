# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines functions that generate workflows for Spin-Orbit calculations.
"""

from fireworks import Firework, Workflow

from matmethods.vasp.firetasks.glue_tasks import PassCalcLocs, CopyVaspOutputs
from matmethods.vasp.firetasks.parse_outputs import VaspToDbTask
from matmethods.vasp.firetasks.run_calc import RunVaspDirect
from matmethods.vasp.firetasks.write_inputs import WriteVaspFromIOSet, ModifyIncar
from matmethods.vasp.input_sets import StructureOptimizationVaspInputSet
from matmethods.utils.utils import soc_warning

__author__ = 'Kiran Mathew'
__credits__ = 'Anubhav jain'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov'


soc_warning()

def get_wf_spinorbit_coupling(structure, magmom, field_directions=[[0,0,1]], vasp_input_set=None,
                              vasp_cmd="vasp", vasp_ncl="vasp_ncl", db_file=None):
    """
    Return Spin-Orbit coupling workflow :

    fw1 : write vasp input set for non-magnetic structural relaxation.
          modify incar to remove magmom settings
          run vasp
          pass run location
          database insertion

    fw2 : Copy files from prev dir
          modify incar for non-magnetic static run, CHGCAR file retained
          run vasp
          pass run location
          database insertion

    soc_fws : list of fireworks(one for each field direction) consisting of the following firetasks:
                 copy files(additional files = CHGCAR) from previous run
                 modify incar
                 run vasp with non-collinear binary(vasp_ncl)
                 pass run location
                 database insertion.

    Args:
        structure (Structure): input structure to be relaxed.
        magmom (list): list of magnetic moment values for each site in the structure.
        field_directions (list): list of magnetic directions for which non-scf vasp soc are to
            be run.
        vasp_input_set (DictVaspInputSet): vasp input set.
        vasp_cmd (string): command to run
        vasp_ncl (string): command to run for non-collinear calculations. Require vasp_ncl binary.
        db_file (string): path to file containing the database credentials.

    Returns:
        Workflow
     """

    task_label = "non-magnetic structure optimization"
    t1 = []
    vasp_input_set = vasp_input_set if vasp_input_set else StructureOptimizationVaspInputSet()
    t1.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))
    t1.append(ModifyIncar(incar_dictmod = {"_unset": {"MAGMOM": ""}}))
    t1.append(RunVaspDirect(vasp_cmd=vasp_cmd))
    t1.append(PassCalcLocs(name=task_label))
    t1.append(VaspToDbTask(db_file=db_file, additional_fields={"task_label": task_label}))
    fw1 = Firework(t1, name="{}-{}".format(structure.composition.reduced_formula, task_label))

    task_label = "non-magnetic static scf"
    t2 = []
    t2.append(CopyVaspOutputs(calc_loc=True))
    t2.append(ModifyIncar(incar_update={"NSW":0, "LCHARG": True}))
    t2.append(RunVaspDirect(vasp_cmd=vasp_cmd))
    t2.append(PassCalcLocs(name=task_label))
    t2.append(VaspToDbTask(db_file=db_file, additional_fields={"task_label": task_label}))
    fw2 = Firework(t2, parents=fw1, name="{}-{}".format(structure.composition.reduced_formula, task_label))

    soc_fws = []
    if len(structure) != len(magmom):
        raise ValueError
    for saxis in field_directions:
        task_label = "non-scf soc " + "".join(str(x) for x in saxis)
        fw_name = "{}-{}".format(structure.composition.reduced_formula, task_label)
        soc_task = []
        additional_files = ["CHGCAR"]
        soc_task.append(CopyVaspOutputs(calc_loc=True, additional_files=additional_files))
        soc_task.append(ModifyIncar(incar_update={"MAGMOM": [[0, 0, m] for m in  magmom],
                                                  "ISYM": -1,
                                                  "LSORBIT": "T",
                                                  "ICHARG": 11,
                                                  "SAXIS": saxis} ))
        soc_task.append(RunVaspDirect(vasp_cmd=vasp_ncl))
        soc_task.append(PassCalcLocs(name=task_label))
        soc_task.append(VaspToDbTask(db_file=db_file, additional_fields={"task_label": task_label}))
        fw = Firework(soc_task, parents=fw2, name=fw_name)
        soc_fws.append(fw)

    return Workflow([fw1, fw2]+soc_fws, name="SOC-"+structure.composition.reduced_formula)
