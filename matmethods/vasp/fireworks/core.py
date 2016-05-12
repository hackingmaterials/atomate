"""
Defines standardized Fireworks that can be chained easily to perform
various sequences of VASP calculations.
"""

from fireworks import Firework
from pymatgen.io.vasp.sets import MPVaspInputSet, MPStaticSet

from matmethods.vasp.firetasks.glue_tasks import CopyVaspOutputs
from matmethods.common.firetasks.glue_tasks import PassCalcLocs
from matmethods.vasp.firetasks.parse_outputs import VaspToDbTask
from matmethods.vasp.firetasks.run_calc import RunVaspDirect
from matmethods.vasp.firetasks.write_inputs import WriteVaspFromIOSet, \
    WriteVaspStaticFromPrev, WriteVaspNSCFFromPrev, \
    WriteVaspDFPTDielectricFromPrev


class OptimizeFW(Firework):
    def __init__(self, structure, name="structure optimization",
                 vasp_input_set=None, vasp_cmd="vasp",
                 override_default_vasp_params=None,
                 db_file=None, parents=None, **kwargs):
        """
        Standard structure optimization Firework.

        Args:
            structure (Structure): Input structure.
            name (str): Name for the Firework.
            vasp_input_set (VaspInputSet): input set to use.
                Defaults to MPVaspInputSet() if None.
            override_default_vasp_params (dict): If this is not None,
                these params are passed to the default vasp_input_set, i.e.,
                MPVaspInputSet. This allows one to easily override some
                settings, e.g., user_incar_settings, etc.
            vasp_cmd (str): Command to run vasp.
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework.
                FW or list of FWS.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        override_default_vasp_params = override_default_vasp_params or {}
        vasp_input_set = vasp_input_set or MPVaspInputSet(
            force_gamma=True, **override_default_vasp_params)

        t = []
        t.append(WriteVaspFromIOSet(structure=structure,
                                    vasp_input_set=vasp_input_set))
        t.append(RunVaspDirect(vasp_cmd=vasp_cmd))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDbTask(db_file=db_file,
                              additional_fields={"task_label": name}))
        super(OptimizeFW, self).__init__(
                t, parents=parents, name="{}-{}".
                format(structure.composition.reduced_formula, name), **kwargs)


class StaticFW(Firework):
    def __init__(self, structure, name="static", vasp_cmd="vasp",
                 copy_vasp_outputs=True, db_file=None, parents=None, **kwargs):
        """
        Standard static calculation Firework.

        Args:
            structure (Structure): Input structure.
            name (str): Name for the Firework.
            vasp_cmd (str): Command to run vasp.
            copy_vasp_outputs (bool): Whether to copy outputs from previous run. Defaults to True.
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        t = []

        if parents:
            if copy_vasp_outputs:
                t.append(
                    CopyVaspOutputs(calc_loc=True, contcar_to_poscar=True))
            t.append(WriteVaspStaticFromPrev(prev_calc_dir='.'))
        else:
            vasp_input_set = MPStaticSet(structure)
            t.append(WriteVaspFromIOSet(structure=structure,
                                        vasp_input_set=vasp_input_set))

        t.append(RunVaspDirect(vasp_cmd=vasp_cmd))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDbTask(db_file=db_file,
                              additional_fields={"task_label": name}))
        super(StaticFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula,
            name), **kwargs)


class NonSCFFW(Firework):
    def __init__(self, structure, name="nscf", mode="uniform", vasp_cmd="vasp",
                 copy_vasp_outputs=True, db_file=None, parents=None, **kwargs):
        """
        Standard NonSCF Calculation Firework supporting both
        uniform and line modes.

        Args:
            structure (Structure): Input structure.
            name (str): Name for the Firework.
            mode (str): uniform or line mode.
            vasp_cmd (str): Command to run vasp.
            copy_vasp_outputs (bool): Whether to copy outputs from previous
                run. Defaults to True.
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework.
                FW or list of FWS.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        t = []
        if copy_vasp_outputs:
            t.append(
                CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"]))
        mode = mode.lower()
        if mode == "uniform":
            t.append(WriteVaspNSCFFromPrev(prev_calc_dir=".", mode="uniform",
                                           reciprocal_density=1000))
        else:
            t.append(WriteVaspNSCFFromPrev(prev_calc_dir=".", mode="line",
                                           reciprocal_density=20))
        t.append(RunVaspDirect(vasp_cmd=vasp_cmd))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDbTask(db_file=db_file, additional_fields={
            "task_label": name + " " + mode},
                              parse_dos=(mode == "uniform"),
                              bandstructure_mode=mode))
        super(NonSCFFW, self).__init__(t, parents=parents, name="%s-%s %s" % (
            structure.composition.reduced_formula,
            name, mode), **kwargs)


class LepsFW(Firework):
    def __init__(self, structure, name="static dielectric", vasp_cmd="vasp",
                 copy_vasp_outputs=True, db_file=None, parents=None, **kwargs):
        """
        Standard static calculation Firework for dielectric constants
        using DFPT.

        Args:
            structure (Structure): Input structure.
            name (str): Name for the Firework.
            vasp_cmd (str): Command to run vasp.
            copy_vasp_outputs (bool): Whether to copy outputs from previous
                run. Defaults to True.
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework.
                FW or list of FWS.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        t = []
        if copy_vasp_outputs:
            t.append(
                CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"],
                                contcar_to_poscar=True))
        t.extend([
            WriteVaspDFPTDielectricFromPrev(prev_calc_dir="."),
            RunVaspDirect(vasp_cmd=vasp_cmd),
            PassCalcLocs(name=name),
            VaspToDbTask(db_file=db_file,
                         additional_fields={"task_label": name})])
        super(LepsFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula,
            name), **kwargs)
