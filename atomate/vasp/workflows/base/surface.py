## for Surface Energy Calculation
from __future__ import division, unicode_literals

"""
#TODO: Write module doc.
       Clean up
"""

__author__ = "Richard Tran"
__version__ = "0.1"
__email__ = "rit634989@gmail.com"
__date__ = "6/24/15"

try:
    # New Py>=3.5 import
    from math import gcd
except ImportError:
    # Deprecated import from Py3.5 onwards.
    from fractions import gcd

import os
import copy

from pymongo import MongoClient

db = MongoClient().data

from custodian.vasp.jobs import VaspJob
from custodian.vasp.handlers import VaspErrorHandler, NonConvergingErrorHandler, \
    UnconvergedErrorHandler, PotimErrorHandler, PositiveEnergyErrorHandler, \
    FrozenJobErrorHandler

from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from atomate.vasp.fireworks.core import OptimizeFW

from pymatgen.core.surface import \
    get_symmetrically_distinct_miller_indices, generate_all_slabs
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.matproj.rest import MPRester
from pymatgen.analysis.structure_analyzer import VoronoiConnectivity
from pymatgen import Element
from pymatgen.io.vasp.sets import MVLSlabSet

from matgendb import QueryEngine

from fireworks.core.firework import Firework, Workflow
from fireworks.core.launchpad import LaunchPad

from matgendb import QueryEngine
from pymongo import MongoClient


class SurfaceDBQueryEngine(QueryEngine):
    def __init__(self, apikey, db_file,
                 tasks_coll="surface_tasks",
                 prop_coll="surface_properties"):
        surf_db_credentials = {'host': dbconfig['host'],
                               'port': dbconfig['port'],
                               'user': dbconfig['user'],
                               'password': dbconfig['password'],
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
        self.db_file = db_file

        super(SurfaceDBQueryEngine, self).__init__(**surf_db_credentials)

    def get_conventional_ucell(formula_id, symprec=1e-3,
                               angle_tolerance=5):
        """
        Gets the conventional unit cell by querying
        materials project for the primitive unit cell

        Args:
            formula_id (string): Materials Project ID
                associated with the slab data entry.
        """

        entries = self.mprester.get_entries(formula_id, inc_structure="Final",
                                            property_data=["e_above_hull"])
        prim_unit_cell = entries[0].structure

        spa = SpacegroupAnalyzer(prim_unit_cell, symprec=symprec,
                                 angle_tolerance=angle_tolerance)
        return spa.get_conventional_standard_structure()


def get_wflow_from_mpid(mpid, qe, vasp_cmd, scratch_dir, db_file,
                        launchpad_dir, max_errors=10,
                        handler_group=[], k_product=50):
    cwd = os.getcwd()
    # Check for conventional unit cell optimization
    conv_ucell_entry = qe.get_entries({"material_id": mpid,
                                       "structure_type": "conventional_unit_cell"},
                                      inc_structure="Final")[0]
    el = conv_ucell_entry.structure[0].species_string

    fws = []
    OptimizeFW_kwargs = {"vasp_cmd": vasp_cmd,
                         "handler_group": handler_group,
                         "scratch_dir": scratch_dir,
                         "job_type": "double_relaxation_run",
                         "half_kpts_first_relax": True,
                         "max_errors": max_errors}
    VaspToDb_kwargs = {"db_file": db_file}
    additional_fields = {"material_id": mpid,
                         "author": os.environ.get("USER")}

    if not conv_ucell_entry:
        ucell = qe.get_conventional_ucell(mpid)
        mvl = MVLSlabSet(ucell, k_product=k_product, bulk=True)
        name = "%s_%s_conventional_k%s" % (el, mpid, k_product)

        OptimizeFW_kwargs["structure"] = ucell
        OptimizeFW_kwargs["vasp_input_set"] = mvl
        OptimizeFW_kwargs["name"] = name
        additional_fields["structure_type"] = "conventional_unit_cell"
        additional_fields["calculation_nme"] = name
        VaspToDb_kwargs["calc_dir"] = os.path.join(cwd, name)
        VaspToDb_kwargs["additional_fields"] = additional_fields
        tasks = [OptimizeFW(**OptimizeFW_kwargs),
                 VaspToDb(**VaspToDb_kwargs)]

    fws.append(tasks)
    wf = Workflow(fws, name='Surface Calculations')
    launchpad = LaunchPad.from_file(os.path.join(os.environ["HOME"],
                                                 launchpad_dir,
                                                 "my_launchpad.yaml"))
    launchpad.add_wf(wf)
