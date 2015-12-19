from fireworks import Firework, Workflow
from matmethods.vasp.firetasks.parse_outputs import VaspToDBTask
from matmethods.vasp.firetasks.run_calc import RunVaspDirect
from matmethods.vasp.firetasks.write_inputs import WriteVaspFromIOSet

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


def get_basic_workflow(structure, vasp_input_set="MPVaspInputSet"):

    write_task = WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set)
    run_task = RunVaspDirect(vasp_cmd="vasp")
    parse_task = VaspToDBTask()

    my_fw = Firework([write_task, run_task, parse_task])

    return Workflow.from_Firework(my_fw)


if __name__ == "__main__":
    pass
