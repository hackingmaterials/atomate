"""
Defines standardized Fireworks that can be chained easily to perform
various sequences of VASP calculations.
"""

from fireworks import Firework
from fireworks import FireTaskBase, explicit_serialize 
from pymatgen.io.vasp.sets import MPVaspInputSet, MPStaticSet, MPSOCSet, MITMDVaspInputSet
from pymatgen.io.vasp.sets import MPRelaxSet, MITMDSet

from matmethods.vasp.firetasks.glue_tasks import CopyVaspOutputs
from matmethods.common.firetasks.glue_tasks import PassCalcLocs
from matmethods.vasp.firetasks.parse_outputs import VaspToDbTask
from matmethods.vasp.firetasks.run_calc import RunVaspCustodian
from matmethods.vasp.firetasks.write_inputs import *
from matmethods.vasp.firetasks.parse_outputs import BandGapCut


class OptimizeFW(Firework):
    def __init__(self, structure, name="structure optimization",
                 vasp_input_set=None, vasp_cmd="vasp",
                 override_default_vasp_params=None, ediffg=None,
                 db_file=None, parents=None, **kwargs):
        """
        Standard structure optimization Firework.

        Args:
            structure (Structure): Input structure.
            name (str): Name for the Firework.
            vasp_input_set (VaspInputSet): input set to use.
                Defaults to MPRelaxSet() if None.
            override_default_vasp_params (dict): If this is not None,
                these params are passed to the default vasp_input_set, i.e.,
                MPRelaxSet. This allows one to easily override some
                settings, e.g., user_incar_settings, etc.
            vasp_cmd (str): Command to run vasp.
            ediffg (float): Shortcut to set ediffg in certain jobs
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework.
                FW or list of FWS.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        override_default_vasp_params = override_default_vasp_params or {}
        vasp_input_set = vasp_input_set or \
                         MPRelaxSet(structure, force_gamma=True,
                                    **override_default_vasp_params)

        t = []
        t.append(WriteVaspFromIOSet(structure=structure,
                                    vasp_input_set=vasp_input_set))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd,
                                  job_type="double_relaxation_run",
                                  max_force_threshold=0.25,
                                  ediffg=ediffg,
                                  auto_npar=">>auto_npar<<"))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDbTask(db_file=db_file,
                              additional_fields={"task_label": name}))
        super(OptimizeFW, self).__init__(
                t, parents=parents, name="{}-{}".
                format(structure.composition.reduced_formula, name), **kwargs)


class StaticFW(Firework):
    def __init__(self, structure, name="static", vasp_input_set=None,
                 vasp_cmd="vasp", copy_vasp_outputs=True, db_file=None,
                 parents=None, calc_loc=True,**kwargs):
        """
        Standard static calculation Firework.

        Args:
            structure (Structure): Input structure.
            name (str): Name for the Firework.
            vasp_input_set (VaspInputSet): input set to use (for jobs w/no parents)
                Defaults to MPStaticSet() if None.
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
                    CopyVaspOutputs(calc_loc=calc_loc, contcar_to_poscar=True))
            t.append(WriteVaspStaticFromPrev(prev_calc_dir='.'))
        else:
            vasp_input_set = MPStaticSet(structure) or vasp_input_set
            t.append(WriteVaspFromIOSet(structure=structure,
                                        vasp_input_set=vasp_input_set))

        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd,
                                  auto_npar=">>auto_npar<<"))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDbTask(db_file=db_file,
                              additional_fields={"task_label": name}))
        super(StaticFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula,
            name), **kwargs)

class HSEBSFW(Firework):
    def __init__(self, structure, parents, name="hse gap", vasp_cmd="vasp",
                 db_file=None, **kwargs):
        """
        For getting a more accurate band gap with HSE - requires previous
        calculation that gives VBM/CBM info. Note that this method is not
        intended for energies, etc. due to sparse k-mesh.

        Args:
            structure (Structure): Input structure.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            name (str): Name for the Firework.
            vasp_cmd (str): Command to run vasp.
            db_file (str): Path to file specifying db credentials.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        t = []
        t.append(CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"]))
        t.append(WriteVaspHSEBSFromPrev(prev_calc_dir='.'))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDbTask(db_file=db_file,
                              additional_fields={"task_label": name}))
        super(HSEBSFW, self).__init__(t, parents=parents, name="{}-{}".format(
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
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd,
                                  auto_npar=">>auto_npar<<"))
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
        if parents:
            if copy_vasp_outputs:
                t.append(
                    CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"],
                                    contcar_to_poscar=True))
                t.append(WriteVaspStaticFromPrev(prev_calc_dir=".",
                                                 lepsilon=True))
        else:
            vasp_input_set = MPStaticSet(structure, lepsilon=True)
            t.append(WriteVaspFromIOSet(structure=structure,
                                        vasp_input_set=vasp_input_set))

        t.extend([
            RunVaspCustodian(vasp_cmd=vasp_cmd),
            PassCalcLocs(name=name),
            VaspToDbTask(db_file=db_file,
                         additional_fields={"task_label": name})])

        super(LepsFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)


class SOCFW(Firework):
    def __init__(self, structure, magmom, name="spinorbit coupling",
                 saxis=(0, 0, 1), vasp_cmd="vasp_ncl",
                 copy_vasp_outputs=True, db_file=None, parents=None, **kwargs):
        """
        Firework for spin orbit coupling calculation.

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

        if parents:
            if copy_vasp_outputs:
                t.append(
                    CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"],
                                    contcar_to_poscar=True))
            t.append(WriteVaspSOCFromPrev(prev_calc_dir=".", magmom=magmom,
                                          saxis=saxis))
        else:
            vasp_input_set = MPSOCSet(structure)
            t.append(WriteVaspFromIOSet(structure=structure,
                                        vasp_input_set=vasp_input_set))
        t.extend([
            RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"),
            PassCalcLocs(name=name),
            VaspToDbTask(db_file=db_file,
                         additional_fields={"task_label": name})])
        super(SOCFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula,
            name), **kwargs)


class TransmuterFW(Firework):
    def __init__(self, structure, transformations, transformation_params=None,
                 vasp_input_set="MPStaticSet", name="structure transmuter", vasp_cmd="vasp",
                 copy_vasp_outputs=True, db_file=None, parents=None, **kwargs):
        """
        Apply the transformations to the input structure, write the input set corresponding
        to the transformed structure and run vasp on them.

        Args:
            structure (Structure): Input structure.
            transformations (list): list of names of transformation classes as defined in
                the modules in pymatgen.transformations. 
                eg:  transformations=['DeformStructureTransformation', 'SupercellTransformation']
            transformation_params (list): list of dicts where each dict specify the input parameters to
                instantiate the transformation class in the transforamtions list.
            vasp_input_set (string): string name for the VASP input set (e.g.,
                "MPStaticSet").
            name (string): Name for the Firework.
            vasp_cmd (string): Command to run vasp.
            copy_vasp_outputs (bool): Whether to copy outputs from previous run. Defaults to True.
            db_file (string): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        t = []

        if parents:
            if copy_vasp_outputs:
                t.append(CopyVaspOutputs(calc_loc=True, contcar_to_poscar=True))
            t.append(WriteTransmutedStructureIOSet(structure=structure, transformations=transformations,
                                                   transformation_params=transformation_params,
                                                   vasp_input_set=vasp_input_set,
                                                   vasp_input_params=kwargs.get("vasp_input_params",{}),
                                                   prev_calc_dir="."))
        else:
            t.append(WriteTransmutedStructureIOSet(structure=structure, transformations=transformations,
                                                   transformation_params=transformation_params,
                                                   vasp_input_set=vasp_input_set,
                                                   vasp_input_params=kwargs.get("vasp_input_params",{})))
        if "vasp_input_params" in kwargs:
            kwargs.pop("vasp_input_params")
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDbTask(db_file=db_file,
                              additional_fields={"task_label": name,
                                                 "transmuter":{"transformations":transformations,
                                                               "transformation_params":transformation_params
                                                              }
                                                }))
        super(TransmuterFW, self).__init__(t, parents=parents,
                                           name="{}-{}".format(structure.composition.reduced_formula, name),
                                           **kwargs)


class MDFW(Firework):
    def __init__(self, structure, start_temp, end_temp, nsteps,
                 name="molecular dynamics run", vasp_input_set=None, vasp_cmd="vasp",
                 override_default_vasp_params=None, wall_time=19200,
                 db_file=None, parents=None, copy_vasp_outputs=True, **kwargs):
        """
        Standard firework for a single MD run.
        Args:
            structure (Structure): Input structure.
            start_temp (float): Start temperature of MD run.
            end_temp (float): End temperature of MD run.
            nsteps (int): Number of MD steps
            name (string): Name for the Firework.
            vasp_input_set (string): string name for the VASP input set (e.g.,
                "MITMDVaspInputSet").
            vasp_cmd (string): Command to run vasp.
            override_default_vasp_params (dict): If this is not None,
                these params are passed to the default vasp_input_set, i.e.,
                MITMDSet. This allows one to easily override some
                settings, e.g., user_incar_settings, etc. Particular to MD,
                one can control time_step and all other settings of the input set.
            wall_time (int): Total wall time in seconds.
            copy_vasp_outputs (bool): Whether to copy outputs from previous run. Defaults to True.
            db_file (string): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """

        override_default_vasp_params = override_default_vasp_params or {}
        vasp_input_set = vasp_input_set or MITMDSet(
            structure, start_temp=start_temp, end_temp=end_temp,
            nsteps=nsteps, **override_default_vasp_params)

        t = []
        if parents:
            if copy_vasp_outputs:
                t.append(CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"],
                                        contcar_to_poscar=True))
        t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, gamma_vasp_cmd=">>gamma_vasp_cmd<<",
                                  handler_group="md", wall_time=wall_time))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDbTask(db_file=db_file,
                              additional_fields={"task_label": name}))
        super(MDFW, self).__init__(
                t, parents=parents, name="{}-{}".
                format(structure.composition.reduced_formula, name), **kwargs)

class LcalcpolFW(Firework):
    def __init__(self, structure, name="static dipole moment", static_name="static", vasp_cmd="vasp",
                 copy_vasp_outputs=False, vasp_input_set=None, db_file=None, parents=None, calc_loc = True,
                 gap_threshold=0.010 **kwargs):
        """
        Standard static calculation Firework for dipole moment. The calculation will not calculate the polarization
        if the band gap of the SCF calculation is metallic (have a band gap less than the gap_threshold).

        The SCF calculation can be provided as a previous run or can be computed within this Firework.

        Args:
            structure (Structure): Input structure.
            name (str): Name for the polarization FireWork.
            static_name (str): Name for the SCF run to be used in PassCalcLoc if copy_vasp_outputs != True.
            vasp_cmd (str): Command to run vasp.
            copy_vasp_outputs (bool): Whether to copy outputs from previous run. Defaults to False.
            vasp_input_set (str): string name for the VASP input set (e.g., "MITMDVaspInputSet").
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            calc_loc (str or True): Name of the previous SCF calculation to be use for CopyVaspOutputs.
                True defaults to parent.
            gap_threshold: Band gap cutoff for determining whether polarization calculation will proceed from
                SCF band gap.
        """
        """
        Standard static calculation Firework for dipole moment. The calculation will not calculate the polarization
        if the band gap of the SCF calculation is metallic (have a band gap less than the gap_threshold).

        The SCF calculation can be provided as a previous run or can be computed within this Firework.

        Args:
            structure (Structure): Input structure.
            name (str): Name for the Firework.
            vasp_cmd (str): Command to run vasp.
            copy_vasp_outputs (bool): Whether to copy outputs from previous
                run. Defaults to True.
            db_file (str): Path to file specifying db credentials.
            parents (Firework): Parents of this particular Firework.
                FW or list of FWS.
            gap_threshold (float): band gap cutoff for determining whether polarization calculation will proceed from
                SCF band gap.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        t = []

        if copy_vasp_outputs:
            t.append(CopyVaspOutputs(calc_loc=calc_loc, additional_files=["CHGCAR","WAVECAR"],contcar_to_poscar=True))

            # Exit Firework if bandgap is less than gap_threshold
            t.append(BandGapCut(gt=gap_threshold,exit_firework=True))
        else:
            # Run a static calculation prior to polarization calculation.
            static = StaticFW(structure, name = static_name, vasp_input_set = vasp_input_set,
                     vasp_cmd = vasp_cmd, copy_vasp_outputs = True, db_file = db_file, parents = parents,
                     calc_loc = cal_loc, ** kwargs)
            t.extend(static.tasks)

            # Exit Firework if bandgap is less than gap_threshold
            t.append(BandGapCut(gt=gap_threshold,exit_firework=True))

            # Create new directory and move to that directory to perform polarization calculation
            polarization_folder = "/polarization"
            if not os.path.exists(os.getcwd() + polarization_folder):
                os.makedirs(os.getcwd() + polarization_folder)
            os.chdir(os.getcwd() + polarization_folder)

            t.append(CopyVaspOutputs(calc_loc=static_name, additional_files=["CHGCAR", "WAVECAR"],
                                     contcar_to_poscar=True))

        t.extend([
            WriteVaspStaticFromPrev(prev_calc_dir=".",other_params={'lcalcpol':True}),
            RunVaspDirect(vasp_cmd=vasp_cmd),
            PassCalcLocs(name=name),
            VaspToDbTask(db_file=db_file,
                         additional_fields={"task_label": name})])

        # Note, Outcar must have read_lcalcpol method for polarization information to be processed.
        # ...assuming VaspDrone will automatically assimilate all properties of the Outcar.

        super(LcalcpolFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula,
            name), **kwargs)
