# coding: utf-8


## for Surface Energy Calculation
from __future__ import division, unicode_literals

__author__ = "Richard Tran"
__version__ = "0.2"
__email__ = "rit001@eng.ucsd.edu"
__date__ = "11/30/17"

try:
    # New Py>=3.5 import
    from math import gcd
except ImportError:
    # Deprecated import from Py3.5 onwards.
    from fractions import gcd

import os
import copy
import json

from pymongo import MongoClient

db = MongoClient().data

from custodian.vasp.jobs import VaspJob
from custodian.vasp.handlers import VaspErrorHandler, NonConvergingErrorHandler, \
    UnconvergedErrorHandler, PotimErrorHandler, PositiveEnergyErrorHandler, \
    FrozenJobErrorHandler

from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from atomate.vasp.fireworks.core import OptimizeFW
from atomate.vasp.firetasks.parse_outputs import VaspToDb

from pymatgen.core.surface import \
    get_symmetrically_distinct_miller_indices, generate_all_slabs
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen import MPRester
from pymatgen.analysis.structure_analyzer import VoronoiConnectivity
from pymatgen import Element
from pymatgen.io.vasp.sets import MVLSlabSet
from pymatgen.io.vasp.outputs import Incar, Outcar

from matgendb import QueryEngine

from fireworks.core.firework import Firework, Workflow, FiretaskBase, FWAction
from fireworks import explicit_serialize
from fireworks.core.launchpad import LaunchPad

from matgendb import QueryEngine
from matgendb.creator import VaspToDbTaskDrone
from pymongo import MongoClient


class SurfacePropertiesWF(object):
    """
    Strings together a bunch of SurfPropFWs (each of these FWs is
    associated with a single structure). Inside SurfPropFWs is a
    bunch of SurfCalcFWs which can represent calculating, inserting
    and post processing of any conv_ucell, ouc, slab, or static slab
    (for work function).
    """

    def __init__(self, apikey, db_file, tasks_coll="surface_tasks",
                 prop_coll="surface_properties", production_mode=False,
                 scratch_dir="", k_product=45, vasp_cmd="vasp", ediffg=-0.02,):

        with open(db_file) as data_file:
            dbconfig = json.load(data_file)

        self.k_product = k_product
        self.vasp_cmd = vasp_cmd
        self.ediffg = ediffg
        self.dbconfig = dbconfig
        self.db_file = db_file
        self.qe = SurfaceDBQueryEngine(dbconfig, apikey)
        self.scratch_dir = scratch_dir

    def wf_from_mpid(self, mpid, polymorph):

        ucell = self.qe.mprester.get_entry_by_material_id(mpid, inc_structure=True,
                                                          conventional_unit_cell=True).structure
        name = "%s_conventional_unit_cell_k%s" % (mpid, self.k_product)

        # bulk = True if calc_type == "oriented_unit_cell" else False
        # get_wf = True if calc_type != "oriented_unit_cell" else False
        mvl = MVLSlabSet(ucell, k_product=self.k_product, bulk=True)

        # This will give us our four basic tasks: WriteVaspFromIOSet,
        # RunVaspCustodian, PassCalcLocs and VaspToDB. We can then
        # modify or remove tasks to suit our needs.
        optimizeFW = OptimizeFW(ucell, name=name, vasp_input_set=mvl,
                                vasp_cmd=self.vasp_cmd, force_gamma=True, parents=None,
                                override_default_vasp_params=None, ediffg=self.ediffg,
                                auto_npar=">>auto_npar<<", job_type="double_relaxation_run")

        tasks = optimizeFW.tasks
        tasks[1] = RunVaspCustodian(vasp_cmd=self.vasp_cmd,
                                    auto_npar=">>auto_npar<<",
                                    job_type="double_relaxation_run",
                                    scratch_dir=self.scratch_dir)

        # Modify list of tasks. This will give us our five tasks:
        # WriteVaspFromIOSet, RunVaspCustodian, PassCalcLocs
        # SurfCalcToDbTask and SurfPropToDbTask
        additional_fields = {"author": os.environ.get("USER"),
                             "structure_type": "conventional_unit_cell",
                             "calculation_name": name}
        # Add mpid as optional so we won't get None
        # when looking for mpid of isolated atoms
        additional_fields["conventional_spacegroup"] = \
            SpacegroupAnalyzer(ucell).get_space_group_symbol()
        additional_fields["polymorph"] = polymorph
        additional_fields["initial_structure"] = ucell
        if mpid:
            additional_fields["material_id"] = mpid

        tasks[3] = VaspToDb(additional_fields=additional_fields, db_file=self.db_file)
        optimizeFW.tasks = tasks

        return Workflow([optimizeFW])


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


class SurfPropToDbTask(FiretaskBase):
    pass

class SurfaceDBQueryEngine(QueryEngine):
    def __init__(self, dbconfig, apikey, tasks_coll="surface_tasks",
                 prop_coll="surface_properties"):

        surf_db_credentials = {'host': dbconfig['host'],
                               'port': dbconfig['port'],
                               'user': dbconfig['admin_user'],
                               'password': dbconfig['admin_password'],
                               'database': dbconfig['database'],
                               'collection': tasks_coll}

        conn = MongoClient(host=surf_db_credentials["host"],
                           port=surf_db_credentials["port"])
        db = conn.get_database(surf_db_credentials["database"])
        db.authenticate(surf_db_credentials["user"],
                        surf_db_credentials["password"])

        surface_properties = db[prop_coll]

        self.surface_properties = surface_properties
        self.mprester = MPRester(apikey)

        super(SurfaceDBQueryEngine, self).__init__(**surf_db_credentials)










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
