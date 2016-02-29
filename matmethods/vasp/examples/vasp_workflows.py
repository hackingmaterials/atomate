from fireworks import Firework, Workflow
from matmethods.vasp.firetasks.glue_tasks import PassVaspLocs, CopyVaspOutputs
from matmethods.vasp.firetasks.parse_outputs import VaspToDBTask
from matmethods.vasp.firetasks.run_calc import RunVaspDirect
from matmethods.vasp.firetasks.write_inputs import WriteVaspFromIOSet, WriteVaspStaticFromPrev
from pymatgen import Lattice, IStructure

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


def get_wf_single_Vasp(structure, vasp_input_set="MPVaspInputSet", vasp_cmd="vasp", db_file=None):

    write_task = WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set)
    run_task = RunVaspDirect(vasp_cmd=vasp_cmd)
    parse_task = VaspToDBTask(db_file=db_file)

    my_fw = Firework([write_task, run_task, parse_task], name="structure optimization")

    return Workflow.from_Firework(my_fw)


def get_wf_double_Vasp(structure, vasp_input_set="MPVaspInputSet", vasp_cmd="vasp", db_file=None):
    t11 = WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set)
    t12 = RunVaspDirect(vasp_cmd=vasp_cmd)
    t13 = PassVaspLocs(name="structure optimization")
    t14 = VaspToDBTask(db_file=db_file, additional_fields={"task_label": "structure optimization"})

    fw1 = Firework([t11, t12, t13, t14], name="structure optimization")

    t21 = CopyVaspOutputs(vasp_loc=True)
    t22 = WriteVaspStaticFromPrev()
    t23 = RunVaspDirect(vasp_cmd=vasp_cmd)
    t24 = PassVaspLocs(name="static")
    t25 = VaspToDBTask(db_file=db_file, additional_fields={"task_label": "static"})

    fw2 = Firework([t21, t22, t23, t24, t25], parents=fw1, name="static")

    return Workflow([fw1, fw2])


if __name__ == "__main__":
    coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
    lattice = Lattice([[3.8401979337, 0.00, 0.00],
                                [1.9200989668, 3.3257101909, 0.00],
                                [0.00, -2.2171384943, 3.1355090603]])
    structure = IStructure(lattice, ["Si"] * 2, coords)

    wf = get_wf_single_Vasp(structure)