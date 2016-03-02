from fireworks import Firework, Workflow
from matmethods.vasp.firetasks.glue_tasks import PassVaspLocs, CopyVaspOutputs
from matmethods.vasp.firetasks.parse_outputs import VaspToDBTask
from matmethods.vasp.firetasks.run_calc import RunVaspDirect
from matmethods.vasp.firetasks.write_inputs import WriteVaspFromIOSet, WriteVaspStaticFromPrev, WriteVaspNSCFFromPrev
from pymatgen import Lattice, IStructure

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


def get_wf_single_Vasp(structure, vasp_input_set="MPVaspInputSet", vasp_cmd="vasp", db_file=None, name="single VASP"):

    write_task = WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set)
    run_task = RunVaspDirect(vasp_cmd=vasp_cmd)
    parse_task = VaspToDBTask(db_file=db_file)

    my_fw = Firework([write_task, run_task, parse_task], name=name)

    return Workflow.from_Firework(my_fw)


def get_wf_bandstructure_Vasp(structure, vasp_input_set="MPVaspInputSet", vasp_cmd="vasp", db_file=None):

    # structure optimization
    t1 = []
    t1.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))
    t1.append(RunVaspDirect(vasp_cmd=vasp_cmd))
    t1.append(PassVaspLocs(name="structure optimization"))
    t1.append(VaspToDBTask(db_file=db_file, additional_fields={"task_label": "structure optimization"}))
    fw1 = Firework(t1, name="structure optimization")

    # static
    t2 = []
    t2.append(CopyVaspOutputs(vasp_loc=True))
    t2.append(WriteVaspStaticFromPrev())
    t2.append(RunVaspDirect(vasp_cmd=vasp_cmd))
    t2.append(PassVaspLocs(name="static"))
    t2.append(VaspToDBTask(db_file=db_file, additional_fields={"task_label": "static"}))
    fw2 = Firework(t2, parents=fw1, name="static")

    # uniform
    t3 = []
    t3.append(CopyVaspOutputs(vasp_loc=True))
    t3.append(WriteVaspNSCFFromPrev(mode="uniform"))
    t3.append(RunVaspDirect(vasp_cmd=vasp_cmd))
    t3.append(PassVaspLocs(name="nscf uniform"))
    t3.append(VaspToDBTask(db_file=db_file, additional_fields={"task_label": "nscf uniform"}, parse_dos=True, bandstructure_mode="uniform"))
    fw3 = Firework(t3, parents=fw2, name="nscf uniform")

    # line mode (run in parallel to uniform)
    t4 = []
    t4.append(CopyVaspOutputs(vasp_loc=True))
    t4.append(WriteVaspNSCFFromPrev(mode="line"))
    t4.append(RunVaspDirect(vasp_cmd=vasp_cmd))
    t4.append(PassVaspLocs(name="nscf line"))
    t4.append(VaspToDBTask(db_file=db_file, additional_fields={"task_label": "nscf line"}, bandstructure_mode="line"))
    fw4 = Firework(t4, parents=fw2, name="nscf line")

    return Workflow([fw1, fw2, fw3, fw4])


if __name__ == "__main__":
    coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
    lattice = Lattice([[3.8401979337, 0.00, 0.00],
                                [1.9200989668, 3.3257101909, 0.00],
                                [0.00, -2.2171384943, 3.1355090603]])
    structure = IStructure(lattice, ["Si"] * 2, coords)

    wf = get_wf_single_Vasp(structure)