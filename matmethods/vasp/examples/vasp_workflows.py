from fireworks import Firework, Workflow
from matmethods.vasp.firetasks.glue_tasks import PassVaspLocs
from matmethods.vasp.firetasks.parse_outputs import VaspToDBTask
from matmethods.vasp.firetasks.run_calc import RunVaspDirect
from matmethods.vasp.firetasks.write_inputs import WriteVaspFromIOSet
from pymatgen import Lattice, IStructure

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


def get_wf_single_Vasp(structure, vasp_input_set="MPVaspInputSet", vasp_cmd="vasp", db_file=None):

    write_task = WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set)
    run_task = RunVaspDirect(vasp_cmd=vasp_cmd)
    parse_task = VaspToDBTask(db_file=db_file)

    my_fw = Firework([write_task, run_task, parse_task])

    return Workflow.from_Firework(my_fw)


def get_wf_double_Vasp(structure, vasp_input_set="MPVaspInputSet", vasp_cmd="vasp", db_file=None):
    t11 = WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set)
    t12 = RunVaspDirect(vasp_cmd=vasp_cmd)
    t13 = PassVaspLocs(name="Structure Optimization")
    t14 = VaspToDBTask(db_file=db_file)

    fw1 = Firework([t11, t12, t13, t14])

    # TODO: t21 should be a STATIC run that reads the previous vasp loc
    # TODO: add a copy task and modify the faker so that it knows that runVasp is now the 3rd task
    t21 = WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set)
    t22 = RunVaspDirect(vasp_cmd=vasp_cmd)
    t23 = PassVaspLocs(name="Static")
    t24 = VaspToDBTask(db_file=db_file)

    fw2 = Firework([t21, t22, t23, t24], parents=fw1)

    return Workflow([fw1, fw2])


if __name__ == "__main__":
    coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
    lattice = Lattice([[3.8401979337, 0.00, 0.00],
                                [1.9200989668, 3.3257101909, 0.00],
                                [0.00, -2.2171384943, 3.1355090603]])
    structure = IStructure(lattice, ["Si"] * 2, coords)

    wf = get_wf_single_Vasp(structure)