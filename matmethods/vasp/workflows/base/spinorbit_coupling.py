# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines functions that generate workflows for Spin-Orbit calculations.
"""

from fireworks import Firework, Workflow

from matmethods.vasp.firetasks.glue_tasks import PassCalcLocs, CopyVaspOutputs
from matmethods.vasp.firetasks.parse_outputs import VaspToDbTask
from matmethods.vasp.firetasks.run_calc import RunVaspDirect
from matmethods.vasp.firetasks.write_inputs import WriteVaspFromIOSet, WriteVaspStaticFromPrev
from matmethods.vasp.input_sets import StructureOptimizationVaspInputSet

__author__ = 'Kiran Mathew'
__credits__ = 'Anubhav jain'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov'


def get_wf_spinorbit_coupling(structure, vasp_input_set=None, vasp_cmd="vasp", db_file=None):
    """
    Return Spin-Orbit coupling workflow consisting of 2 fireworks :

    Firework 1 : write vasp input set for structural relaxation. Retain CHGCAR and WAVECAR files
                 run vasp,
                 pass run location,
                 database insertion.

    Firework 2 : copy files(additional files = CHGCAR, WAVECAR) from previous run,
                 write vasp input set for static run with incar settings override,
                 run vasp with noncllinear binary(vasp_ncl)
                 pass run location
                 database insertion.

    Args:
        structure (Structure): input structure to be relaxed.
        vasp_input_set (DictVaspInputSet): vasp input set.
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.

    Returns:
        Workflow
     """

    task_label = "structure optimization"
    t1 = []
    # need WAVECAR and CHGCAR for the next step
    config_dict_override = {"INCAR": {"LWAVE": "T", "LCHARG": "T"}}
    vasp_input_set = vasp_input_set if vasp_input_set \
        else StructureOptimizationVaspInputSet(config_dict_override=config_dict_override)

    t1.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set,
                                 config_dict_override=config_dict_override))
    t1.append(RunVaspDirect(vasp_cmd=vasp_cmd))
    t1.append(PassCalcLocs(name=task_label))
    t1.append(VaspToDbTask(db_file=db_file,
                           additional_fields={"task_label": task_label}))
    fw1 = Firework(t1, name="{}-{}".format(structure.composition.reduced_formula, task_label))

    task_label = "non-scf soc"
    t2 = []
    config_dict_override = {"INCAR": {"ISYM": -1,
                                      "LSORBIT": "T",
                                      "LNONCOLLINEAR": "T",
                                      "ICHARG": 11,
                                      "SAXIS": [0, 0, 1]}}
    preserve_magmom = True
    preserve_old_incar = True
    additional_files = ["CHGCAR", "WAVECAR"]
    t2.append(CopyVaspOutputs(calc_loc=True, additional_files=additional_files))
    t2.append(WriteVaspStaticFromPrev(config_dict_override=config_dict_override,
                                      preserve_magmom=preserve_magmom,
                                      preserve_old_incar=preserve_old_incar))
    t2.append(RunVaspDirect(vasp_cmd="srun vasp_ncl"))
    t2.append(PassCalcLocs(name=task_label))
    t2.append(VaspToDbTask(db_file=db_file, additional_fields={"task_label": task_label}))
    fw2 = Firework(t2, parents=fw1, name="{}-{}".format(structure.composition.reduced_formula,
                                                        task_label))

    return Workflow([fw1, fw2], name=structure.composition.reduced_formula)
