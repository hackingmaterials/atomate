# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import json
import os
import shutil
import unittest
import zlib

import gridfs
from pymongo import MongoClient, DESCENDING

from fireworks import LaunchPad, FWorker
from fireworks.core.rocket_launcher import rapidfire

from matmethods.vasp.vasp_powerups import use_custodian, add_namefile, use_fake_vasp, add_trackers
from matmethods.vasp.workflows.base.core import get_wf
from matmethods.vasp.workflows.presets.core import wf_raman_spectra

from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.util.testing import PymatgenTest


__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "reference_files", "db_connections")
ref_dir = os.path.join(module_dir, "test_files")

DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...


class TestVaspWorkflows(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not os.environ.get("VASP_PSP_DIR"):
            os.environ["VASP_PSP_DIR"] = os.path.join(module_dir, "..", "..", "tests", "reference_files")
            print(
                'Note: This system is not set up to run VASP jobs. '
                'Please set your VASP_PSP_DIR environment variable.')

        cls.struct_si = PymatgenTest.get_structure("Si")
        cls.scratch_dir = os.path.join(module_dir, "scratch")
        cls.raman_config = {"modes": [0, 1], "step_size": 0.005,
                            "vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<"}
        cls.wf = wf_raman_spectra(cls.struct_si, cls.raman_config)

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
        reference_dir = os.path.abspath(os.path.join(ref_dir, "raman_wf"))
        si_ref_dirs = {"structure optimization": os.path.join(reference_dir, "1"),
                       "phonon static dielectric": os.path.join(reference_dir, "2"),
                       "raman_0_-0.005 static dielectric": os.path.join(reference_dir, "6"),
                       "raman_0_0.005 static dielectric": os.path.join(reference_dir, "5"),
                       "raman_1_-0.005 static dielectric": os.path.join(reference_dir, "4"),
                       "raman_1_0.005 static dielectric": os.path.join(reference_dir, "3"),
                       "raman analysis": os.path.join(reference_dir, "7"),}
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
        if mode not in ["structure optimization", "phonon static dielectric",
                        "raman_0_0.005 static dielectric", "raman analysis"]:
            raise ValueError("Invalid mode!")

        if mode not in ["raman analysis"]:
            self.assertEqual(d["formula_pretty"], "Si")
            self.assertEqual(d["formula_anonymous"], "A")
            self.assertEqual(d["nelements"], 1)
            self.assertEqual(d["state"], "successful")
            self.assertAlmostEqual(d["calcs_reversed"][0]["output"]["structure"]["lattice"]["a"], 3.867, 2)

        if mode in ["structure optimization"]:
            self.assertAlmostEqual(d["output"]["energy"], -10.850, 2)
            self.assertAlmostEqual(d["output"]["energy_per_atom"], -5.425, 2)
        elif mode in ["phonon static dielectric"]:
            # TODO
            pass
        elif mode in ["raman_0_0.005 static dielectric"]:
            # TODO
            pass
        elif mode in ["raman analysis"]:
            # TODO
            pass

    def test_wf(self):
        self.wf = self._simulate_vasprun(self.wf)

        self.assertEqual(len(self.wf.fws), len(self.raman_config["modes"]) * 2 + 3)

        self.lp.add_wf(self.wf)
        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))

        d = self._get_task_collection().find_one({"task_label": "structure optimization"})
        self._check_run(d, mode="structure optimization")

        # make sure the static run ran OK
        d = self._get_task_collection().find_one({"task_label": "phonon static dielectric"})
        self._check_run(d, mode="phonon static dielectric")

        # make sure the uniform run ran OK
        d = self._get_task_collection().find_one({"task_label": "raman_0_0.005 static dielectric"})
        self._check_run(d, mode="raman_0_0.005 static dielectric")

        # make sure the uniform run ran OK
        d = self._get_task_collection().find_one({"task_label": "raman analysis"})
        self._check_run(d, mode="raman analysis")


if __name__ == "__main__":
    unittest.main()
