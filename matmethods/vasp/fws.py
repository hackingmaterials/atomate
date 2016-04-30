from pymatgen.io.vasp.sets import MPVaspInputSet, MPStaticSet

from fireworks import Firework

from matmethods.vasp.firetasks.glue_tasks import PassCalcLocs, CopyVaspOutputs
from matmethods.vasp.firetasks.parse_outputs import VaspToDbTask
from matmethods.vasp.firetasks.run_calc import RunVaspDirect
from matmethods.vasp.firetasks.write_inputs import WriteVaspFromIOSet, \
    WriteVaspStaticFromPrev, WriteVaspNSCFFromPrev, \
    WriteVaspDFPTDielectricFromPrev


class OptimizeFW(Firework):
    def __init__(self, structure, name="structure optimization",
                 vasp_input_set=None, vasp_cmd="vasp",
                 db_file=None, parents=None, **kwargs):
        vasp_input_set = vasp_input_set if vasp_input_set else MPVaspInputSet(
            force_gamma=True)

        t = []
        t.append(WriteVaspFromIOSet(structure=structure,
                                    vasp_input_set=vasp_input_set))
        t.append(RunVaspDirect(vasp_cmd=vasp_cmd))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDbTask(db_file=db_file,
                              additional_fields={"task_label": name}))
        super(OptimizeFW, self).__init__(t, parents=parents,
                                         name="{}-{}".format(
                                             structure.composition.reduced_formula,
                                             name), **kwargs)


class StaticFW(Firework):
    def __init__(self, structure, name="static", vasp_cmd="vasp",
                 copy_vasp_outputs=True,
                 db_file=None, parents=None, **kwargs):

        t = []

        if parents:
            if copy_vasp_outputs:
                t.append(CopyVaspOutputs(calc_loc=True, contcar_to_poscar=True))
            t.append(WriteVaspStaticFromPrev())
        else:
            vasp_input_set = MPStaticSet()
            t.append(WriteVaspFromIOSet(structure=structure,
                                        vasp_input_set=vasp_input_set))

        t.append(RunVaspDirect(vasp_cmd=vasp_cmd))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDbTask(db_file=db_file,
                              additional_fields={"task_label": name}))
        super(StaticFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula,
            name), **kwargs)


class NonSCFUniformFW(Firework):
    def __init__(self, structure, name="nscf uniform", vasp_cmd="vasp",
                 copy_vasp_outputs=True,
                 db_file=None, parents=None, **kwargs):
        t = []
        if copy_vasp_outputs:
            t.append(
                CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"]))
        t.append(WriteVaspNSCFFromPrev(mode="uniform"))
        t.append(RunVaspDirect(vasp_cmd=vasp_cmd))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDbTask(db_file=db_file,
                              additional_fields={"task_label": name},
                              parse_dos=True, bandstructure_mode="uniform"))
        super(NonSCFUniformFW, self).__init__(t, parents=parents,
                                              name="{}-{}".format(
                                                  structure.composition.reduced_formula,
                                                  name), **kwargs)


class NonSCFLineFW(Firework):
    def __init__(self, structure, name="nscf line", vasp_cmd="vasp",
                 copy_vasp_outputs=True,
                 db_file=None, parents=None, **kwargs):
        t = []
        if copy_vasp_outputs:
            t.append(
                CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"]))
        t.append(WriteVaspNSCFFromPrev(mode="line"))
        t.append(RunVaspDirect(vasp_cmd=vasp_cmd))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDbTask(db_file=db_file,
                              additional_fields={"task_label": name},
                              bandstructure_mode="line"))
        super(NonSCFLineFW, self).__init__(t, parents=parents,
                                           name="{}-{}".format(
                                               structure.composition.reduced_formula,
                                               name), **kwargs)


class LepsFW(Firework):
    def __init__(self, structure, name="static dielectric", vasp_cmd="vasp",
                 copy_vasp_outputs=True,
                 db_file=None, parents=None, **kwargs):
        t = []
        if copy_vasp_outputs:
            t.append(
                CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"], contcar_to_poscar=True))
        t.extend([
             WriteVaspDFPTDielectricFromPrev(),
             RunVaspDirect(vasp_cmd=vasp_cmd),
             PassCalcLocs(name=name),
             VaspToDbTask(db_file=db_file,
                          additional_fields={"task_label": name})])
        super(LepsFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula,
            name), **kwargs)
