from fireworks import Firework

from matmethods.vasp.firetasks.glue_tasks import PassCalcLocs, CopyVaspOutputs
from matmethods.vasp.firetasks.parse_outputs import VaspToDbTask
from matmethods.vasp.firetasks.run_calc import RunVaspDirect
from matmethods.vasp.firetasks.write_inputs import WriteVaspFromIOSet, \
    WriteVaspStaticFromPrev, WriteVaspNSCFFromPrev, WriteVaspDFPTDielectricFromPrev
from matmethods.vasp.input_sets import StructureOptimizationVaspInputSet
from matmethods.vasp.vasp_powerups import decorate_write_name
from matmethods.vasp.workflows.base.single_vasp import get_wf_single


class OptimizeFW(Firework):

    def __init__(self, structure, name="structure optimization", vasp_input_set=None, vasp_cmd="vasp",
                 db_file=None, parents=None, **kwargs):
        vasp_input_set = vasp_input_set if vasp_input_set else StructureOptimizationVaspInputSet()

        t1 = []
        t1.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))
        t1.append(RunVaspDirect(vasp_cmd=vasp_cmd))
        t1.append(PassCalcLocs(name=name))
        t1.append(VaspToDbTask(db_file=db_file,
                               additional_fields={"task_label": name}))
        super(OptimizeFW, self).__init__(t1, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)


class StaticFW(Firework):

    def __init__(self, structure, name="static", vasp_cmd="vasp", copy_vasp_outputs=True,
                 db_file=None, parents=None, **kwargs):

        t2 = []
        if copy_vasp_outputs:
            t2.append(CopyVaspOutputs(calc_loc=True))
        t2.append(WriteVaspStaticFromPrev())
        t2.append(RunVaspDirect(vasp_cmd=vasp_cmd))
        t2.append(PassCalcLocs(name=name))
        t2.append(VaspToDbTask(db_file=db_file, additional_fields={"task_label": name}))
        super(StaticFW, self).__init__(t2, parents=parents, name="{}-{}".format(structure.composition.reduced_formula,
                                           name), **kwargs)


class NonSCFUniformFW(Firework):

    def __init__(self, structure, name="nscf uniform", vasp_cmd="vasp", copy_vasp_outputs=True,
                 db_file=None, parents=None, **kwargs):
        t3 = []
        if copy_vasp_outputs:
            t3.append(CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"]))
        t3.append(WriteVaspNSCFFromPrev(mode="uniform"))
        t3.append(RunVaspDirect(vasp_cmd=vasp_cmd))
        t3.append(PassCalcLocs(name=name))
        t3.append(VaspToDbTask(db_file=db_file,
                               additional_fields={"task_label": name},
                               parse_dos=True, bandstructure_mode="uniform"))
        super(NonSCFUniformFW, self).__init__(t3, parents=parents, name="{}-{}".format(structure.composition.reduced_formula,
                                           name), **kwargs)


class NonSCFLineFW(Firework):

    def __init__(self, structure, name="nscf line", vasp_cmd="vasp", copy_vasp_outputs=True,
                 db_file=None, parents=None, **kwargs):
        t4 = []
        if copy_vasp_outputs:
            t4.append(CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"]))
        t4.append(WriteVaspNSCFFromPrev(mode="line"))
        t4.append(RunVaspDirect(vasp_cmd=vasp_cmd))
        t4.append(PassCalcLocs(name=name))
        t4.append(VaspToDbTask(db_file=db_file,
                               additional_fields={"task_label": name},
                               bandstructure_mode="line"))
        super(NonSCFLineFW, self).__init__(t4, parents=parents, name="{}-{}".format(structure.composition.reduced_formula,
                                           name), **kwargs)

class LepsFW(Firework):

    def __init__(self, structure, name="static dielectric", vasp_cmd="vasp", copy_vasp_outputs=True,
                 db_file=None, parents=None, **kwargs):
        t5 = [CopyVaspOutputs(calc_loc=True, contcar_to_poscar=True),
              WriteVaspDFPTDielectricFromPrev(),
              RunVaspDirect(vasp_cmd=vasp_cmd),
              PassCalcLocs(name=name),
              VaspToDbTask(db_file=db_file,
                           additional_fields={"task_label": name})]
        super(LepsFW, self).__init__(t5, parents=parents, name="{}-{}".format(structure.composition.reduced_formula,
                                           name), **kwargs)
