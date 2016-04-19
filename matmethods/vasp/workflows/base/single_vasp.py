from fireworks import Firework, Workflow
from matmethods.vasp.firetasks.parse_outputs import VaspToDbTask, ToDbTask
from matmethods.vasp.firetasks.run_calc import RunVaspDirect
from matmethods.vasp.firetasks.write_inputs import WriteVaspFromIOSet
from matmethods.vasp.input_sets import StructureOptimizationVaspInputSet
from matmethods.vasp.vasp_powerups import use_custodian
from matmethods.vasp.drones import VaspToDbTaskDrone

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


def get_wf_single_Vasp(structure, vasp_input_set=None, vasp_cmd="vasp",
                       db_file=None, task_label="single VASP",
                       custodian_powerup=False):
    """
    Return vasp workflow consisting of a single firework. The firework is made
    of  3 firetasks: write inputset for structural relaxation, run vasp and
    database insertion

    Args:
        structure (Structure): input structure to be relaxed.
        vasp_input_set (DictVaspInputSet): vasp input set.
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.
        task_label (str): workflow name.
        custodian_powerup (bool): If set the vasp run will be wrapped in
            custodian.

    Returns:
        Workflow
    """
    vasp_input_set = vasp_input_set if vasp_input_set else StructureOptimizationVaspInputSet()

    write_task = WriteVaspFromIOSet(structure=structure,
                                    vasp_input_set=vasp_input_set)
    run_task = RunVaspDirect(vasp_cmd=vasp_cmd)
    parse_task = ToDbTask(db_file=db_file, drone=VaspToDbTaskDrone())

    my_fw = Firework([write_task, run_task, parse_task],
                     name="{}-{}".format(structure.composition.reduced_formula,
                                         task_label))
    my_wf = Workflow.from_Firework(my_fw)

    if custodian_powerup:
        my_wf = use_custodian(my_wf)

    return my_wf
