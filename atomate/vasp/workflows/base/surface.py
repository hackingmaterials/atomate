# coding: utf-8


## for Surface Energy Calculation
from __future__ import division, unicode_literals

__author__ = "Richard Tran"
__version__ = "0.2"
__email__ = "rit001@eng.ucsd.edu"
__date__ = "11/30/17"

import json

from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from atomate.vasp.firetasks.parse_outputs import VaspToDb
from atomate.common.firetasks.glue_tasks import CreateFolder
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet

from pymatgen.core.surface import generate_all_slabs, Structure, SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.sets import MVLSlabSet

from fireworks.core.firework import Firework, Workflow, FiretaskBase, FWAction
from fireworks import explicit_serialize


class SurfacePropertiesWF(object):
    """
    Strings together a bunch of SurfPropFWs (each of these FWs is
    associated with a single structure). Inside SurfPropFWs is a
    bunch of SurfCalcFWs which can represent calculating, inserting
    and post processing of any conv_ucell, ouc, slab, or static slab
    (for work function).
    """

    def __init__(self, db_file, tasks_coll="surface_tasks",
                 prop_coll="surface_properties", production_mode=False,
                 scratch_dir="", k_product=50, vasp_cmd="vasp"):

        with open(db_file) as data_file:
            dbconfig = json.load(data_file)

        self.k_product = k_product
        self.vasp_cmd = vasp_cmd
        self.dbconfig = dbconfig
        self.db_file = db_file
        self.scratch_dir = scratch_dir
        self.production_mode = production_mode
        self.tasks_coll = tasks_coll
        self.prop_coll = prop_coll

    def wf_from_mpid(self, structure, mmi, mpid=None):

        return Workflow([ConvUcellFW(structure, mmi, self.tasks_coll, self.prop_coll,
                                     self.production_mode, self.scratch_dir,
                                     self.k_product, self.db_file, self.vasp_cmd, mpid=mpid)])

class ConvUcellFW(Firework):
    def __init__(self, ucell, mmi, tasks_coll, prop_coll,
                 production_mode, scratch_dir, k_product,
                 db_file, vasp_cmd, mpid="--", **kwargs):
        """
        Customized FW similar to OptimizeFW.

        Args:
            tasks_coll (str): Name of the collection to store the raw data
            prop_coll (str): Name of the collection to store the post-processed
                data (to be queried from the MAPI)
            production_mode (bool): Whether or not to run addition fws, tasks,
                operations meant for high-throughput production, set to False for
                most users who just run a workflow for non high-throughput purposes.
            scratch_dir: (str) - if specified, uses this directory as the root
                scratch dir. Supports env_chk.
            k_product (int): Default to 50, kpoint number * length for a & b
                directions, also for c direction in bulk calculations.
            db_file (str): FULL path to file containing the database credentials.
                Supports env_chk.
            vasp_cmd (str): Command to run vasp.
            structure (Structure): Input structure. Need to input either the
                structure or an mpid
            mpid (str): Unique Materials Project structure ID of the unit cell.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """

        name = "%s_%s_conventional_unit_cell_k%s" % \
               (ucell[0].species_string, mpid, k_product)
        mvl = MVLSlabSet(ucell, k_product=k_product, bulk=True)

        # This will give us our four basic tasks: WriteVaspFromIOSet,
        # RunVaspCustodian, PassCalcLocs and VaspToDB. We can then
        # modify or remove tasks to suit our needs.
        tasks = []
        tasks.append(CreateFolder(folder_name=name, change_dir=True))
        tasks.append(WriteVaspFromIOSet(structure=ucell, vasp_input_set=mvl))
        tasks.append(RunVaspCustodian(vasp_cmd=vasp_cmd,
                                      auto_npar=">>auto_npar<<",
                                      job_type="double_relaxation_run",
                                      scratch_dir=scratch_dir))

        # Add additional fields to distinguish this calculation
        additional_fields = {"structure_type": "conventional_unit_cell",
                             "calculation_name": name}
        additional_fields["conventional_spacegroup"] = \
            SpacegroupAnalyzer(ucell).get_space_group_symbol()
        additional_fields["initial_structure"] = ucell
        additional_fields["material_id"] = mpid
        tasks.append(VaspToDb(additional_fields=additional_fields, db_file=db_file))
        tasks.append(FacetFWsGeneratorTask(structure_type="conventional_unit_cell", mmi=mmi,
                                           vasp_cmd=vasp_cmd, mpid=mpid, db_file=db_file,
                                           scratch_dir=scratch_dir, k_product=k_product,
                                           tasks_coll=tasks_coll, prop_coll=prop_coll,
                                           production_mode=production_mode))

        super(ConvUcellFW, self).__init__(tasks, name=name, **kwargs)

class OUCFW(Firework):
    def __init__(self, ouc, tasks_coll, prop_coll, miller_index, scale_factor,
                 production_mode, scratch_dir, k_product,
                 db_file, vasp_cmd, mpid="--", **kwargs):
        """
        Customized FW similar to OptimizeFW.

        Args:
            tasks_coll (str): Name of the collection to store the raw data
            prop_coll (str): Name of the collection to store the post-processed
                data (to be queried from the MAPI)
            production_mode (bool): Whether or not to run addition fws, tasks,
                operations meant for high-throughput production, set to False for
                most users who just run a workflow for non high-throughput purposes.
            scratch_dir: (str) - if specified, uses this directory as the root
                scratch dir. Supports env_chk.
            k_product (int): Default to 50, kpoint number * length for a & b
                directions, also for c direction in bulk calculations.
            db_file (str): FULL path to file containing the database credentials.
                Supports env_chk.
            vasp_cmd (str): Command to run vasp.
            structure (Structure): Input structure. Need to input either the
                structure or an mpid
            mpid (str): Unique Materials Project structure ID of the unit cell.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """

        name = "%s_%s_bulk_k%s_%s%s%s" % \
               (ouc[0].species_string, mpid, k_product,
                miller_index[0], miller_index[1], miller_index[2])
        mvl = MVLSlabSet(ouc, k_product=k_product, bulk=True)

        # This will give us our four basic tasks: WriteVaspFromIOSet,
        # RunVaspCustodian, PassCalcLocs and VaspToDB. We can then
        # modify or remove tasks to suit our needs.
        tasks = []
        tasks.append(CreateFolder(folder_name=name, change_dir=True))
        tasks.append(WriteVaspFromIOSet(structure=ouc, vasp_input_set=mvl))
        tasks.append(RunVaspCustodian(vasp_cmd=vasp_cmd,
                                      auto_npar=">>auto_npar<<",
                                      job_type="double_relaxation_run",
                                      scratch_dir=scratch_dir))

        # Add additional fields to distinguish this calculation
        additional_fields = {"structure_type": "oriented_unit_cell",
                             "calculation_name": name}
        additional_fields["conventional_spacegroup"] = \
            SpacegroupAnalyzer(ouc).get_space_group_symbol()
        additional_fields["initial_structure"] = ouc
        additional_fields["material_id"] = mpid
        additional_fields["miller_index"] = miller_index
        additional_fields["scale_factor"] = scale_factor
        tasks.append(VaspToDb(additional_fields=additional_fields, db_file=db_file))

        tasks.append(FacetFWsGeneratorTask(structure_type="oriented_unit_cell",
                                           vasp_cmd=vasp_cmd, mpid=mpid, db_file=db_file,
                                           scratch_dir=scratch_dir, k_product=k_product,
                                           tasks_coll=tasks_coll, prop_coll=prop_coll,
                                           production_mode=production_mode,
                                           miller_index=miller_index))

        super(OUCFW, self).__init__(tasks, name=name, **kwargs)

class SlabFW(Firework):
    def __init__(self, slab, tasks_coll, prop_coll, slab_gen_params,
                 miller_index, scale_factor, ouc, shift, ssize, vsize,
                 production_mode, scratch_dir, k_product, db_file,
                 vasp_cmd, reconstruction=None, mpid="--", **kwargs):
        """
        Customized FW similar to OptimizeFW.

        Args:
            tasks_coll (str): Name of the collection to store the raw data
            prop_coll (str): Name of the collection to store the post-processed
                data (to be queried from the MAPI)
            production_mode (bool): Whether or not to run addition fws, tasks,
                operations meant for high-throughput production, set to False for
                most users who just run a workflow for non high-throughput purposes.
            scratch_dir: (str) - if specified, uses this directory as the root
                scratch dir. Supports env_chk.
            k_product (int): Default to 50, kpoint number * length for a & b
                directions, also for c direction in bulk calculations.
            db_file (str): FULL path to file containing the database credentials.
                Supports env_chk.
            vasp_cmd (str): Command to run vasp.
            structure (Structure): Input structure. Need to input either the
                structure or an mpid
            mpid (str): Unique Materials Project structure ID of the unit cell.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """

        if not reconstruction:
            name = "%s_%s_slab_k%s_s%sv%s_%s%s%s_shift%s" % \
                   (slab[0].species_string, mpid, k_product, ssize, vsize,
                    miller_index[0], miller_index[1], miller_index[2], shift)
        else:
            name = "%s_%s_slab_k%s_s%sv%s_%s" % (slab[0].species_string, mpid, \
                                                 k_product, ssize, vsize, reconstruction)

        mvl = MVLSlabSet(slab, k_product=k_product, bulk=False, get_wf=True)

        # This will give us our four basic tasks: WriteVaspFromIOSet,
        # RunVaspCustodian, PassCalcLocs and VaspToDB. We can then
        # modify or remove tasks to suit our needs.
        tasks = []
        tasks.append(CreateFolder(folder_name=name, change_dir=True))
        tasks.append(WriteVaspFromIOSet(structure=slab, vasp_input_set=mvl))
        tasks.append(RunVaspCustodian(vasp_cmd=vasp_cmd,
                                      auto_npar=">>auto_npar<<",
                                      job_type="double_relaxation_run",
                                      scratch_dir=scratch_dir))

        # Add additional fields to distinguish this calculation
        additional_fields = {"structure_type": "slab_cell",
                             "calculation_name": name}
        additional_fields["conventional_spacegroup"] = \
            SpacegroupAnalyzer(ouc).get_space_group_symbol()
        additional_fields["initial_structure"] = slab
        additional_fields["material_id"] = mpid
        additional_fields["miller_index"] = miller_index
        additional_fields["scale_factor"] = scale_factor
        additional_fields["oriented_unit_cell"] = ouc
        additional_fields["slab_size"] = ssize
        additional_fields["vac_size"] = vsize
        additional_fields["shift"] = shift
        additional_fields["reconstruction"] = reconstruction
        # additional_fields["local_potential_along_c"]
        # additional_fields["efermi"]
        tasks.append(VaspToDb(additional_fields=additional_fields, db_file=db_file))

        super(SlabFW, self).__init__(tasks, name=name, **kwargs)


@explicit_serialize
class FacetFWsGeneratorTask(FiretaskBase):
    """
    Task for generating FWs for oriented unit cell calculations or slab
    calculations. If the initial structure is an oriented unit cell,
    generate slab fws, if its a conventional unit cell, generate oriented
    unit cell and reconstruction calculations.
    """

    required_params = ['structure_type', "scratch_dir", "k_product",
                       "db_file", "tasks_coll", "prop_coll", "vasp_cmd"]
    optional_params = ["slab_gen_params", "production_mode", "mpid", "mmi", "miller_index"]

    def run_task(self, fw_spec):

        if self.get("slab_gen_params", None):
            slab_gen_params = self.get("slab_gen_params", None)
        else:
            slab_gen_params = {"min_slab_size": 10, "min_vacuum_size": 10,
                               "symmetrize": True, "center_slab": True}
            slab_gen_params["max_normal_search"] = self.get("mmi") if \
                self.get("mmi") else max(self.get("miller_index"))
            if self.get("structure_type") == "conventional_unit_cell":
                slab_gen_params["include_reconstructions"] = True

        FWs = []
        if self.get("structure_type") == "conventional_unit_cell":
            # Then we create a set of FWs for oriented_unit_cell
            # calculations and reconstructed slabs

            slab_gen_params["structure"] = Structure.from_file("CONTCAR.relax2.gz")

            all_slabs = generate_all_slabs(max_index=self.get("mmi"),
                                           **slab_gen_params)
            miller_list = []
            for slab in all_slabs:
                if slab.reconstruction:
                    FWs.append(self.get_slab_fw(slab, slab_gen_params))
                else:
                    if tuple(slab.miller_index) in miller_list:
                        continue
                    else:
                        # build a oriented unit cell fw
                        FWs.append(self.get_ouc_fw(slab))

        elif self.get("structure_type") == "oriented_unit_cell":
            slab_gen_params["initial_structure"] = Structure.from_file("CONTCAR.relax2.gz")
            slab_gen_params["miller_index"] = [0,0,1]
            symmetrize = slab_gen_params["symmetrize"]
            del slab_gen_params["symmetrize"]
            slabgen = SlabGenerator(**slab_gen_params)
            for slab in slabgen.get_slabs(symmetrize=symmetrize):
                slab.miller_index = self.get("miller_index")
                FWs.append(self.get_slab_fw(slab, slab_gen_params))

        return FWAction(additions=FWs)

    def get_ouc_fw(self, slab):

        return OUCFW(slab.oriented_unit_cell, self.get("tasks_coll"),
                     self.get("prop_coll"), slab.miller_index,
                     slab.scale_factor, self.get("production_mode", True),
                     self.get("scratch_dir"), self.get("k_product"),
                     self.get("db_file"), self.get("vasp_cmd"),
                     mpid=self.get("mpid", "--"))

    def get_slab_fw(self, slab, slab_gen_params):

        return SlabFW(slab, self.get("tasks_coll"),
                      self.get("prop_coll"), slab_gen_params,
                      slab.miller_index, slab.scale_factor,
                      slab.oriented_unit_cell, slab.shift,
                      slab_gen_params["min_slab_size"],
                      slab_gen_params["min_vacuum_size"],
                      self.get("production_mode", True),
                      self.get("scratch_dir"), self.get("k_product"),
                      self.get("db_file"), self.get("vasp_cmd"),
                      reconstruction=slab.reconstruction,
                      mpid=self.get("mpid", "--"))

# @explicit_serialize
# class SurfPropToDbTask(FiretaskBase):
#     """
#     Analyzes the stress/strain data of an elastic workflow to produce
#     an elastic tensor and various other quantities.
#
#     Required params:
#         structure (Structure): structure to use for symmetrization,
#             input structure.  If an optimization was used, will
#             look for relaxed structure in calc locs
#
#     Optional params:
#         db_file (str): path to file containing the database credentials.
#             Supports env_chk. Default: write data to JSON file.
#         order (int): order of fit to perform
#         fw_spec_field (str): if set, will update the task doc with the contents
#             of this key in the fw_spec.
#         fitting_method (str): if set, will use one of the specified
#             fitting methods from pymatgen.  Supported methods are
#             "independent", "pseudoinverse", and "finite_difference."
#             Note that order 3 and higher required finite difference
#             fitting, and will override.
#     """
#
#     required_params = ['structure']
#     optional_params = ['db_file', 'order', 'fw_spec_field', 'fitting_method']
#
#     def run_task(self, fw_spec):
#
#         # Add additional fields to distinguish this calculation
#         additional_fields = {"structure_type": "conventional_unit_cell",
#                              "calculation_name": name}
#         additional_fields["conventional_spacegroup"] = \
#             SpacegroupAnalyzer(ucell).get_space_group_symbol()
#         additional_fields["initial_structure"] = ucell
#         # Add mpid as optional so we won't get None
#         # when looking for mpid of isolated atoms
#         if mpid:
#             additional_fields["material_id"] = mpid
#
#         return FWAction()


# @explicit_serialize
# class SurfCalcToDbTask(FiretaskBase):
#     """
#         Inserts a single vasp calculation in a folder into a DB and
#         useful information pertaining to slabs and oriented unit cells.
#     """
#
#     required_params = ["vaspdbinsert_parameters", "struct_type", "polymorph"]
#     optional_params = ["surface_area", "shift", "debug", "diatomic", "mpid",
#                        "miller_index", "vsize", "ssize", "isolated_atom"]
#
#     def run_task(self, fw_spec):
#         """
#             Required Parameters:
#                 host (str): See SurfaceWorkflowManager in surface_wf.py
#                 port (int): See SurfaceWorkflowManager in surface_wf.py
#                 user (str): See SurfaceWorkflowManager in surface_wf.py
#                 password (str): See SurfaceWorkflowManager in surface_wf.py
#                 database (str): See SurfaceWorkflowManager in surface_wf.py
#                 collection (str): See SurfaceWorkflowManager in surface_wf.py
#                 mpid (str): The Materials Project ID associated with the
#                     initial structure used to build the slab from
#                 struct_type (str): either oriented_unit_cell or slab_cell
#                 miller_index (list): Miller Index of the oriented
#                     unit cell or slab
#                 loc (str path): Location of the outputs of
#                     the vasp calculations
#                 cwd (str): Current working directory
#                 conventional_spacegroup (str): The spacegroup of the structure
#                     asscociated with the MPID input
#                 polymorph (str): The rank of the  polymorph of the structure
#                     associated with the MPID input, 0 being the ground state
#                     polymorph.
#             Optional Parameters:
#                 surface_area (float): surface area of the slab, obtained
#                     from slab object before relaxation
#                 shift (float): A shift value in Angstrom that determines how
#                     much a slab should be shifted. For determining number of
#                     terminations, obtained from slab object before relaxation
#                 vsize (float): Size of vacuum layer of slab in Angstroms,
#                     obtained from slab object before relaxation
#                 ssize (float): Size of slab layer of slab in Angstroms,
#                     obtained from slab object before relaxation
#                 isolated_atom (str): Specie of the structure used to
#                     calculate the energy of an isolated atom (for cohesive
#                     energy calculations)
#         """
#
#         # Get all the optional/required parameters
#         # dec = MontyDecoder()
#         struct_type = self.get("struct_type")
#         shift = self.get("shift", None)
#         vsize = self.get("vsize", None)
#         ssize = self.get("ssize", None)
#         miller_index = self.get("miller_index")
#         mpid = self.get("mpid", None)
#         polymorph = self.get("polymorph")
#         vaspdbinsert_parameters = self.get("vaspdbinsert_parameters")
#
#         warnings = []
#         # Addtional info relating to slabs
#         additional_fields = {"author": os.environ.get("USER"),
#                              "structure_type": struct_type,
#                              "final_incar": Incar.from_file("./INCAR.relax2.gz"),
#                              "final_magnetization": Outcar("./OUTCAR.relax2.gz").magnetization,
#                              "calculation_name": name,
#                              "warnings": warnings}
#
#         # Add mpid as optional so we won't get None
#         # when looking for mpid of isolated atoms
#         additional_fields["miller_index"] = miller_index
#         additional_fields["surface_area"] = surface_area
#         additional_fields["shift"] = shift
#         additional_fields["vac_size"] = vsize
#         additional_fields["slab_size"] = ssize
#         additional_fields["material_id"] = mpid
#         additional_fields["conventional_spacegroup"] = spacegroup
#         additional_fields["polymorph"] = polymorph
#
#         if mpid:
#             additional_fields["material_id"] = mpid
#
#         drone = VaspToDbTaskDrone(use_full_uri=False,
#                                   additional_fields=additional_fields,
#                                   **vaspdbinsert_parameters)
#         drone.assimilate(calc_locs)













        # if self.calc_type == "conventional_unit_cell":
        #     name = "-%s_conventional_unit_cell_k%s" %(self.mpid, self.k_product)
        # elif self.calc_type == "oriented_unit_cell":
        #     "-%s_bulk_k%s_%s%s%s" % (self.mpid, self.k_product, self.hkl[0],
        #                              self.hkl[1], self.hkl[2])
        # elif self.calc_type == "slab_cell":
        #     if reconstruction_name:
        #         return "-%s_slab_k%s_s%sv%s_%s" %(self.mpid, self.k_product,
        #                                           self.min_slab_size,
        #                                           self.min_vac_size,
        #                                           self.reconstruction_name)
        #     else:
        #         return "-%s_slab_k%s_s%sv%s_%s%s%s_shift%s" % (self.mpid, self.k_product,
        #                                                        self.min_slab_size,
        #                                                        self.min_vac_size,
        #                                                        self.hkl[0], self.hkl[1],
        #                                                        self.hkl[2], self.shift)







































# def get_fw_from_ucell(ucell, vasp_cmd, scratch_dir, db_file,
#                       max_errors=10, handler_group=[], k_product=50,
#                       mpid="--", inc_conv_ucell=False):
#
#     # Get the pretty formula for naming
#     comp = ucell.composition.reduced_composition.as_dict()
#     pretty_formula = ""
#     for el, stoich in comp.items():
#         pretty_formula += str(el)
#         if stoich != 1:
#             pretty_formula += str(int(stoich))
#
#     # Set up Firetask parameters
#     OptimizeFW_kwargs = {"vasp_cmd": vasp_cmd,
#                          "handler_group": handler_group,
#                          "scratch_dir": scratch_dir,
#                          "job_type": "double_relaxation_run",
#                          "half_kpts_first_relax": True,
#                          "max_errors": max_errors}
#     VaspToDb_kwargs = {"db_file": db_file}
#     additional_fields = {"pretty_formula": pretty_formula,
#                          "material_id": mpid,
#                          "author": os.environ.get("USER")}
#
#     cwd = os.getcwd()
#     fws = []
#     if inc_conv_ucell:
#         mvl = MVLSlabSet(ucell, k_product=k_product, bulk=True)
#         name = "%s_%s_conventional_k%s" % (pretty_formula, mpid, k_product)
#
#         # Edit Firetask parameters for particular Firework
#         OptimizeFW_kwargs["structure"] = ucell
#         OptimizeFW_kwargs["vasp_input_set"] = mvl
#         OptimizeFW_kwargs["name"] = name
#         additional_fields["structure_type"] = "conventional_unit_cell"
#         additional_fields["calculation_nme"] = name
#         VaspToDb_kwargs["calc_dir"] = os.path.join(cwd, name)
#         VaspToDb_kwargs["additional_fields"] = additional_fields
#
#         tasks = [OptimizeFW(**OptimizeFW_kwargs),
#                  VaspToDb(**VaspToDb_kwargs)]
#         fws.append(Firework(tasks, name="%s_%s" %(pretty_formula, mpid)))
#
#     return fws
#
# def get_wflow_from_mpid(mpid, qe, vasp_cmd, scratch_dir, db_file,
#                         max_errors=10, handler_group=[], k_product=50):
#
#     # Check for conventional unit cell optimization
#     conv_ucell_entry = qe.get_entries({"material_id": mpid,
#                                        "structure_type": "conventional_unit_cell"},
#                                       inc_structure="Final")
#
#     inc_conv_ucell = True if not conv_ucell_entry else False
#     ucell = conv_ucell_entry.structure if conv_ucell_entry \
#         else qe.get_conventional_ucell(mpid)
#
#     fws = get_tasks_from_ucell(ucell, qe, vasp_cmd, scratch_dir,
#                                db_file, max_errors=max_errors,
#                                handler_group=handler_group, k_product=k_product,
#                                mpid=mpid, inc_conv_ucell=inc_conv_ucell)
#
#     return Workflow(fws, name='Surface Calculations')
        # launchpad = LaunchPad.from_file(os.path.join(os.environ["HOME"],
        #                                              launchpad_dir,
        #                                              "my_launchpad.yaml"))
        # launchpad.add_wf(wf)
