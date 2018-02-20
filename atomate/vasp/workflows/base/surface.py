# coding: utf-8


## for Surface Energy Calculation
from __future__ import division, unicode_literals

__author__ = "Richard Tran"
__version__ = "0.2"
__email__ = "rit001@eng.ucsd.edu"
__date__ = "11/30/17"

import json
import os

from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from atomate.common.firetasks.glue_tasks import CreateFolder, RenameFile
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet
from atomate.utils.utils import env_chk
from atomate.vasp.database import VaspCalcDb
from atomate.common.firetasks.glue_tasks import get_calc_loc
from atomate.utils.utils import get_logger
from atomate.vasp.drones import VaspDrone
from atomate.vasp.firetasks.parse_outputs import VaspToDb

from pymatgen.core.surface import generate_all_slabs, Structure, SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.sets import MVLSlabSet
from pymatgen.analysis.surface_analysis import SlabEntry, SurfaceEnergyPlotter, ComputedStructureEntry
from pymatgen import MPRester

from fireworks.core.firework import Firework, Workflow, FiretaskBase, FWAction
from fireworks import explicit_serialize

logger = get_logger(__name__)


class SurfacePropertiesWF(object):
    """
    Strings together a bunch of SurfPropFWs (each of these FWs is
    associated with a single structure). Inside SurfPropFWs is a
    bunch of SurfCalcFWs which can represent calculating, inserting
    and post processing of any conv_ucell, ouc, slab, or static slab
    (for work function).
    """

    def __init__(self, db_file, tasks_coll="surface_tasks", cwd=os.getcwd(),
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
        self.cwd = cwd

    def wf_from_mpid(self, structure, mmi, mpid=None):

        return Workflow([ConvUcellFW(structure, mmi, self.tasks_coll, self.prop_coll,
                                     self.production_mode, self.scratch_dir, self.k_product,
                                     self.db_file, self.vasp_cmd, cwd=self.cwd, mpid=mpid)])

class ConvUcellFW(Firework):
    def __init__(self, ucell, mmi, tasks_coll, prop_coll,
                 production_mode, scratch_dir, k_product,
                 db_file, vasp_cmd, cwd=os.getcwd(), mpid="--", **kwargs):
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
        tasks.append(CreateFolder(folder_name=os.path.join(cwd, name),
                                  change_dir=True, relative_path=True))
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
        additional_fields["initial_structure"] = ucell.as_dict()
        additional_fields["material_id"] = mpid
        tasks.append(VaspToDb(additional_fields=additional_fields, db_file=db_file))
        # tasks.append(SurfPropToDbTask(structure_type="conventional_unit_cell",
        #                               db_file=db_file, additional_fields=additional_fields))
        tasks.append(FacetFWsGeneratorTask(structure_type="conventional_unit_cell", mmi=mmi,
                                           vasp_cmd=vasp_cmd, mpid=mpid, db_file=db_file,
                                           scratch_dir=scratch_dir, k_product=k_product,
                                           tasks_coll=tasks_coll, prop_coll=prop_coll,
                                           production_mode=production_mode, cwd=cwd))

        super(ConvUcellFW, self).__init__(tasks, name=name, **kwargs)

class OUCFW(Firework):
    def __init__(self, ouc, tasks_coll, prop_coll, miller_index,
                 scale_factor, production_mode, scratch_dir, k_product,
                 db_file, vasp_cmd, cwd=os.getcwd(), mpid="--", **kwargs):
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
        tasks.append(CreateFolder(folder_name=os.path.join(cwd, name),
                                  change_dir=True, relative_path=True))
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
        additional_fields["initial_structure"] = ouc.as_dict()
        additional_fields["material_id"] = mpid
        additional_fields["miller_index"] = miller_index
        additional_fields["scale_factor"] = scale_factor
        tasks.append(VaspToDb(additional_fields=additional_fields, db_file=db_file))
        # tasks.append(SurfPropToDbTask(structure_type="oriented_unit_cell",
        #                               db_file=db_file, additional_fields=additional_fields))
        tasks.append(FacetFWsGeneratorTask(structure_type="oriented_unit_cell",
                                           vasp_cmd=vasp_cmd, mpid=mpid, db_file=db_file,
                                           scratch_dir=scratch_dir, k_product=k_product,
                                           tasks_coll=tasks_coll, prop_coll=prop_coll,
                                           production_mode=production_mode,
                                           miller_index=miller_index, cwd=cwd))

        super(OUCFW, self).__init__(tasks, name=name, **kwargs)

class SlabFW(Firework):
    def __init__(self, slab, tasks_coll, prop_coll, slab_gen_params,
                 miller_index, scale_factor, ouc, shift, ssize, vsize,
                 production_mode, scratch_dir, k_product, db_file,
                 vasp_cmd, reconstruction=None, cwd=os.getcwd(),
                 mpid="--", **kwargs):
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
        tasks.append(CreateFolder(folder_name=os.path.join(cwd, name),
                                  change_dir=True, relative_path=True))
        tasks.append(WriteVaspFromIOSet(structure=slab, vasp_input_set=mvl))
        tasks.append(RunVaspCustodian(vasp_cmd=vasp_cmd,
                                      auto_npar=">>auto_npar<<",
                                      job_type="double_relaxation_run",
                                      scratch_dir=scratch_dir))

        # Add additional fields to distinguish this calculation
        additional_fields = {"structure_type": "slab_cell",
                             "calculation_name": name}
        sg = SpacegroupAnalyzer(ouc)
        additional_fields["conventional_spacegroup"] = \
            {"symbol": sg.get_space_group_symbol(), "number": sg.get_space_group_number()}
        additional_fields["initial_structure"] = slab.as_dict()
        additional_fields["material_id"] = mpid
        additional_fields["miller_index"] = miller_index
        additional_fields["scale_factor"] = scale_factor
        additional_fields["oriented_unit_cell"] = ouc.as_dict()
        additional_fields["slab_size"] = ssize
        additional_fields["vac_size"] = vsize
        additional_fields["shift"] = shift
        additional_fields["reconstruction"] = reconstruction
        # additional_fields["local_potential_along_c"]
        # additional_fields["efermi"]
        tasks.append(RenameFile(file="LOCPOT.gz", new_name="LOCPOT.relax2.gz"))
        tasks.append(VaspToDb(additional_fields=additional_fields, db_file=db_file))
        # tasks.append(SurfPropToDbTask(structure_type="slab_cell", prop_coll=prop_coll,
        #                               additional_fields=additional_fields,
        #                               db_file=db_file))
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
    optional_params = ["slab_gen_params", "production_mode", "mpid", "mmi", "miller_index", "cwd"]

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
                      cwd=self.get("cwd", os.getcwd()),
                      mpid=self.get("mpid", "--"))

EV_PER_ANG2_TO_JOULES_PER_M2 = 16.0217656

@explicit_serialize
class SurfPropToDbTask(FiretaskBase):
    """
    Analyzes the stress/strain data of an elastic workflow to produce
    an elastic tensor and various other quantities.

    Required params:
        structure (Structure): structure to use for symmetrization,
            input structure.  If an optimization was used, will
            look for relaxed structure in calc locs

    Optional params:
        db_file (str): path to file containing the database credentials.
            Supports env_chk. Default: write data to JSON file.
        order (int): order of fit to perform
        fw_spec_field (str): if set, will update the task doc with the contents
            of this key in the fw_spec.
        fitting_method (str): if set, will use one of the specified
            fitting methods from pymatgen.  Supported methods are
            "independent", "pseudoinverse", and "finite_difference."
            Note that order 3 and higher required finite difference
            fitting, and will override.
    """
    required_params = ["structure_type"]
    optional_params = ["prop_coll", "additional_fields", "db_file"]

    def run_task(self, fw_spec):

        # get the database connection
        self.db_file = env_chk(self.get('db_file'), fw_spec)
        with open(self.db_file) as db_file:
            self.dbconfig = json.load(db_file)

        self.mmdb = VaspCalcDb.from_db_file(self.db_file, admin=True)
        self.mprester = MPRester()
        self.surftasks = self.mmdb.db[self.dbconfig["collection"]]

        # Insert the raw calculation
        self.task_doc = self.parse_raw_calcs()
        self.mpid = self.task_doc["material_id"]
        self.mpentry = self.mprester.get_entry_by_material_id(self.task_doc["material_id"],
                                                              property_data=["e_above_hull"])

        # If the prop_coll is given, store post processed surface energies in
        # that collection, if not, store it in default collection from db_file
        if self.get("prop_coll", None):
            self.dbconfig["collection"] = self.get("prop_coll", None)
        self.surfpropdb = VaspCalcDb(**self.dbconfig)
        self.surfprops = self.get_surface_properties_entry()

        # # If this is a slab calculation, we should have all
        # # the data needed to get the surface properties now
        # if self.get("structure_type") == "slab_cell":
        #
        #     logger.info("Retrieving surface properties")
        #
        #     propdoc = self.get_facet_properties()
        #
        #     # Get the facet specific properties. Note there can be more than
        #     # one item for a Miller index when reconstruction is a possibility
        #     surfaces, facetprops = [], []
        #     for surface in self.surfprops["surfaces"]:
        #         if tuple(surface["miller_index"]) == tuple(self.task_doc["miller_index"]):
        #             facetprops.append(surface)
        #         else:
        #             surfaces.append(surface)
        #
        #     # Three possibilities:
        #     # 1) there are no facet entries as of yet
        #     # 2) only one entry exists (its unreconstructed)
        #     # 3) two entries exists (its reconstructed)
        #
        #     # Do entries for this specific facet exist already?
        #     if not facetprops:
        #         # If not, then just append the current facet entry
        #         surfaces.append(propdoc)
        #     elif len(facetprops) > 1:
        #         # That means one of these entries must be reconstructed
        #         resurf = [surf for surf in facetprops if surf['is_reconstructed']][0]
        #         surf = [surf for surf in facetprops if not surf['is_reconstructed']][0]
        #         if not propdoc["is_reconstructed"] and \
        #                 all([se < surf["surface_energy"],
        #                      se < resurf["surface_energy"]]):
        #             # then replace both entries, reconstruction won't happen and
        #             # the current unreconstructed termination is too unstable
        #             surfaces.append(propdoc)
        #         elif not propdoc["is_reconstructed"] \
        #                 and propdoc["surface_energy"] < surf["surface_energy"]:
        #             # It is only more stable than the unreconstructed
        #             # entry, replace that entry only
        #             surfaces.append(propdoc)
        #             surfaces.append(resurf)
        #         elif propdoc["is_reconstructed"] and \
        #                         propdoc["surface_energy"] < resurf["surface_energy"]:
        #             # it is a more stable reconstruction, append it along with the unreconstructed
        #             surfaces.append(propdoc)
        #             surfaces.append(surf)
        #         else:
        #             # do nothing
        #             surfaces.extend(facetprops)
        #     else:
        #         # That means no reconstruction as of yet
        #         if propdoc["surface_energy"] < facetprops[0]["surface_energy"]:
        #             # replace the current doc
        #             surfaces.append(propdoc)
        #         else:
        #             # do not do anything
        #             surfaces.extend(facetprops)
        #
        #     # Get properties relative to the Wulff shape (overall surfaces)
        #     self.surfprops["surfaces"] = surfaces
        #     self.surfprops = self.get_relative_properties(self.surfprops)
        #     self.surfpropdb.update_one({"material_id": self.mpid},
        #                                {"$set": self.surfprops})
        #
        #     logger.info("Finished parsing surface properties with task_id: {}".format(t_id))

        if self.get("defuse_unsuccessful", True):
            defuse_children = (self.task_doc["state"] != "successful")
        else:
            defuse_children = False

        return FWAction(stored_data={"task_id": self.task_doc.get("task_id", None)},
                        defuse_children=defuse_children)

    def parse_raw_calcs(self):

        # This just acts as the VaspToDb FireTask

        # get the directory that contains the VASP dir to parse
        calc_dir = os.getcwd()

        # parse the VASP directory
        logger.info("PARSING DIRECTORY: {}".format(calc_dir))

        # get the locpot for work function calculations
        if self.get("structure_type") == "slab_cell":
            os.rename("LOCPOT.gz", "LOCPOT.relax2.gz")

        drone = VaspDrone(additional_fields=self.get("additional_fields"),
                          parse_dos=self.get("parse_dos", False), compress_dos=1,
                          bandstructure_mode=self.get("bandstructure_mode",
                                                      False), compress_bs=1)

        # assimilate (i.e., parse)
        task_doc = drone.assimilate(calc_dir)

        # Check for additional keys to set based on the fw_spec
        if self.get("fw_spec_field"):
            task_doc.update(fw_spec[self.get("fw_spec_field")])

        # db insertion or taskdoc dump
        if not self.db_file:
            with open("task.json", "w") as f:
                f.write(json.dumps(task_doc, default=DATETIME_HANDLER))
        else:
            t_id = self.mmdb.insert_task(task_doc, parse_dos=self.get("parse_dos", False),
                                         parse_bs=bool(self.get("bandstructure_mode", False)))

            logger.info("Finished parsing raw tasks with task_id: {}".format(t_id))

        return task_doc

    def get_facet_properties(self):

        # Get surface properties for a specific facet
        facet_props = {}

        # Get the ouc and slab entries
        ouc_task = self.surftasks.find_one({"material_id": self.mpid,
                                            "structure_type": "oriented_unit_cell",
                                            "miller_index": self.task_doc["miller_index"]})
        slab = Structure.from_dict(self.task_doc["structure"])
        slab_entry = SlabEntry(slab, self.task_doc["energy"], self.task_doc["miller_index"])
        ouc = Structure.from_dict(ouc_task["structure"])
        ouc_entry = ComputedStructureEntry(ouc, ouc_task["energy"])
        se = slab_entry.surface_energy(ouc_entry)
        facet_props["surface_energy_EV_PER_ANG2"] = se
        facet_props["surface_energy"] = se*EV_PER_ANG2_TO_JOULES_PER_M2
        facet_props["miller_index"] = self.task_doc["miller_index"]




        return

    def get_surface_properties_entry(self):

        surfprops = self.surfpropdb.find_one({"material_id": self.mpid})

        if surfprops:
            # Check if a material entry even exists for this system
            return surfprops
        else:
            # Then we create the entry.
            surfprops = {}
            surfprops["surfaces"] = []
            surfprops["e_above_hull"] = self.mpentry.data["e_above_hull"]
            entries = mprester.get_entries(self.surftasks["elements"][0] ,
                                           property_data=["e_above_hull", "material_id"])
            sorted_entries = sorted(entries, key=lambda entry: entry.data["e_above_hull"])
            surfprops["weighted_surface_energy_EV_PER_ANG2"] = None
            surfprops["weighted_surface_energy"] = None
            surfprops["pretty_formula"] = self.surftasks["elements"][0]
            surfprops["material_id"] = self.task_doc["material_id"]
            surfprops["polymorph"] = [i for i, entry in enumerate(sorted_entries)
                                           if entry.data["material_id"] == self.mpid][0]
            surfprops["spacegroup"] = self.task_doc["conventional_spacegroup"]
            surfprops["surface_anisotropy"] = None
            surfprops["shape_factor"] = None
            surfprops["weighted_work_function"] = None

            # Add this new entry to the collection
            self.surfpropdb.insert(surfprops)
            return surfprops


    def get_wulff_shape(self):

        return

    def get_relative_properties(self):

        return