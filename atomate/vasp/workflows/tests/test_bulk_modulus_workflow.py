# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import json
import os
import shutil
import unittest

import numpy as np
from monty.json import MontyEncoder

from pymongo import MongoClient

from fireworks import LaunchPad, FWorker, Workflow
from fireworks.core.rocket_launcher import rapidfire

from atomate.vasp.powerups import use_fake_vasp
from atomate.vasp.workflows.presets.core import wf_bulk_modulus

from pymatgen import SETTINGS, Structure
from pymatgen.util.testing import PymatgenTest
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "reference_files", "db_connections")
reference_dir = os.path.abspath(os.path.join(module_dir, "test_files", "bulk_modulus_wf"))

DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...


class TestBulkModulusWorkflow(unittest.TestCase):
    """This tests can be run in two modes: 1. first all inputs and outputs of all deformations are present
    in which case """
    @classmethod
    def setUpClass(cls):
        if not SETTINGS.get("PMG_VASP_PSP_DIR"):
            SETTINGS["PMG_VASP_PSP_DIR"] = os.path.join(module_dir, "..", "..", "tests", "reference_files")
            print('This system is not set up to run VASP jobs. '
                  'Please set PMG_VASP_PSP_DIR variable in your ~/.pmgrc.yaml file.')

        cls.struct_si = PymatgenTest.get_structure("Si")
        cls.scratch_dir = os.path.join(module_dir, "scratch")
        cls.elastic_config = {"vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<"}
        cls.wf = wf_bulk_modulus(cls.struct_si, strain_max=0.05, nsteps=6)
        cls.pre_run_deformations = [1, 2, 3, 5]



    def setUp(self):
        if os.path.exists(self.scratch_dir):
            shutil.rmtree(self.scratch_dir)
        os.makedirs(self.scratch_dir)
        os.chdir(self.scratch_dir)
        try:
            self.lp = LaunchPad.from_file(os.path.join(db_dir, "my_launchpad.yaml"))
            self.lp.reset("", require_password=False)
        except:
            raise unittest.SkipTest(
                'Cannot connect to MongoDB! Is the database server running? '
                'Are the credentials correct?')



    def tearDown(self):
        if not DEBUG_MODE:
            shutil.rmtree(self.scratch_dir)
            self.lp.reset("", require_password=False)
            db = self._get_task_database()
            for coll in db.collection_names():
                if coll != "system.indexes":
                    db[coll].drop()



    def _simulate_vasprun(self, wf):
        si_ref_dirs = {"bulk_modulus deformation 0": os.path.join(reference_dir, str(2)),
                       "bulk_modulus deformation 4": os.path.join(reference_dir, str(6))}

        # for the first time use_fake_vasp is called to create the json files.
        for i in self.pre_run_deformations:
            if os.path.exists(os.path.join(reference_dir, str(i), "inputs")):
                si_ref_dirs["bulk_modulus deformation {}".format(i)] = os.path.join(reference_dir, str(i+2))


        si_ref_dirs["structure optimization"] =  os.path.join(reference_dir, "1")
        return use_fake_vasp(wf, si_ref_dirs, params_to_check=["ENCUT"])



    def _get_task_database(self):
        with open(os.path.join(db_dir, "db.json")) as f:
            creds = json.loads(f.read())
            conn = MongoClient(creds["host"], creds["port"])
            db = conn[creds["database"]]
            if "admin_user" in creds:
                db.authenticate(creds["admin_user"], creds["admin_password"])
            return db



    def _get_task_collection(self, coll_name=None):
        with open(os.path.join(db_dir, "db.json")) as f:
            creds = json.loads(f.read())
            db = self._get_task_database()
            coll_name = coll_name or creds["collection"]
            return db[coll_name]



    def _check_run(self, d, mode):
        if mode not in ["structure optimization", "bulk_modulus deformation 0",
                        "bulk_modulus deformation 4", "fit equation of state"]:
            raise ValueError("Invalid mode!")

        if mode not in ["fit equation of state"]:
            self.assertEqual(d["formula_pretty"], "Si")
            self.assertEqual(d["formula_anonymous"], "A")
            self.assertEqual(d["nelements"], 1)
            self.assertEqual(d["state"], "successful")

        #TODO: the following number must be corrected, they are just copied from test_elastic_workflow.py
        if mode in ["structure optimization"]:
            self.assertAlmostEqual(d["calcs_reversed"][0]["output"]["structure"]["lattice"]["a"], 3.866, 3)
            self.assertAlmostEqual(d["output"]["energy_per_atom"], -5.432, 3)

        elif mode in ["bulk_modulus deformation 0"]:
            stress = d["calcs_reversed"][0]["output"]["ionic_steps"][-1]["stress"]
            np.testing.assert_allclose(stress, np.diag([189.19, 189.19, 189.19]), atol=1e-2)

        elif mode in ["bulk_modulus deformation 4"]:
            stress = d["calcs_reversed"][0]["output"]["ionic_steps"][-1]["stress"]
            np.testing.assert_allclose(stress, np.diag([-65.56, -65.56, -65.56]), atol=1e-2)
        elif mode in ["fit equation of state"]:
            self.assertAlmostEqual(d["bulk_modulus"], 88.90, places=2)
            self.assertEqual(len(d["all_task_ids"]), 7)
            self.assertEqual(len(d["energies"]), 6)
            self.assertEqual(len(d["volumes"]), 6)
            s = SpacegroupAnalyzer(Structure.from_dict(d["structure"])).get_conventional_standard_structure()
            self.assertAlmostEqual(s.lattice.c, 5.468, places=3)



    def test_wf(self):
        self.wf = self._simulate_vasprun(self.wf)
        self.assertEqual(len(self.wf.fws), 8)

        # check vasp parameters for ionic relaxation
        defo_vis = [fw.spec["_tasks"][2]['vasp_input_set'] for fw in self.wf.fws if "deform" in fw.name]
        assert all([vis['user_incar_settings']['NSW'] == 99 for vis in defo_vis])
        assert all([vis['user_incar_settings']['IBRION'] == 2 for vis in defo_vis])

        self.lp.add_wf(self.wf)
        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))

        # either write the task docs that are not tested to files (only once) or read them to avoid too many files
        task_from_json = False
        for i in self.pre_run_deformations:
            if not os.path.exists(os.path.join(reference_dir, str(i + 2), "inputs")):
                with open(os.path.join(reference_dir, str(i + 2), "deformation_{}.json".format(i))) as fp:
                    d = json.load(fp)
                    new_task = self.lp.fireworks.find_one({"name": {"$regex": "bulk_modulus deformation {}".format(i)}})

                    # make the tag in the old task loaded from .json compatible with the new task
                    d["task_label"] = new_task["name"]

                    # to avoid duplicated task_ids:
                    d["task_id"] = -d["task_id"]
                    db = self._get_task_database()
                    db["tasks"].insert(d)
                task_from_json = True

                # bypass the pre-run deformations if the inputs don't exist; otherwise, write such tasks to file
                self.lp.fireworks.update_one({"name": {"$regex": "bulk_modulus deformation {}".format(i)}},
                                             {'$set': {"state": "COMPLETED"}})

            else:
                # this step needs to be run  once: when deformation*.json is present, remove the inputs/outputs folders
                d = self._get_task_collection().find_one(
                    {"task_label": {"$regex": "bulk_modulus deformation {}".format(i)}})
                with open(os.path.join(reference_dir, str(i + 2), "deformation_{}.json".format(i)), 'w') as fp:
                    json.dump(d, fp, sort_keys=True, indent=4, ensure_ascii=False, cls=MontyEncoder)

        if task_from_json:
            self.lp.fireworks.update_one({"name": {"$regex": "fit equation of state"}}, {'$set': {"state": "READY"}})
            rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))


        # check relaxation
        d = self._get_task_collection().find_one({"task_label": {"$regex": "structure optimization"}})
        self._check_run(d, mode="structure optimization")

        # check two of the deformation calculations
        d = self._get_task_collection().find_one({"task_label": {"$regex": "bulk_modulus deformation 0"}})
        self._check_run(d, mode="bulk_modulus deformation 0")

        d = self._get_task_collection().find_one({"task_label": {"$regex":"bulk_modulus deformation 4"}})
        self._check_run(d, mode="bulk_modulus deformation 4")

        # check the final results
        d = self._get_task_collection(coll_name="eos").find_one()
        self._check_run(d, mode="fit equation of state")