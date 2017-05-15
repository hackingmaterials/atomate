# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import json
import os
import shutil
import unittest

import numpy as np

from pymongo import MongoClient

from fireworks import LaunchPad, FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.vasp.powerups import use_fake_vasp
from atomate.vasp.workflows.base.adsorption import get_wf_adsorption, get_wf_adsorption_from_slab

from pymatgen import SETTINGS, Structure, Molecule, Lattice
from pymatgen.util.testing import PymatgenTest
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

__author__ = 'Kiran Mathew, Joseph Montoya'
__email__ = 'montoyjh@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
ref_dir = os.path.join(module_dir, "..", "..", "test_files")

DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...


class TestAdsorptionWorkflow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not SETTINGS.get("PMG_VASP_PSP_DIR"):
            SETTINGS["PMG_VASP_PSP_DIR"] = os.path.join(module_dir, "..", "..", "tests", "..", "..", "test_files")
            print('This system is not set up to run VASP jobs. '
                  'Please set PMG_VASP_PSP_DIR variable in your ~/.pmgrc.yaml file.')

        cls.struct_ir = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3.875728), ["Ir"], [[0, 0, 0]])
        cls.scratch_dir = os.path.join(module_dir, "scratch")
        cls.ads_config = {"100": [Molecule("H", [[0, 0, 0]])]}
        cls.wf_1 = get_wf_adsorption(cls.struct_ir, cls.ads_config, 
                                     db_file=os.path.join(db_dir, "db.json"))

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
            os.chdir(module_dir)

    def _simulate_vasprun(self, wf):
        reference_dir = os.path.abspath(os.path.join(ref_dir, "adsorbate_wf"))
        ir_ref_dirs = {"Ir-structure optimization": os.path.join(reference_dir, "1"),
                       "Ir-Ir_100 slab optimization": os.path.join(reference_dir, "2"),
                       "Ir-H1-Ir_100 adsorbate optimization 0": os.path.join(reference_dir, "3"),
                       "Ir-H1-Ir_100 adsorbate optimization 1": os.path.join(reference_dir, "4"),
                       "Ir-H1-Ir_100 adsorbate optimization 2": os.path.join(reference_dir, "5")}
        return use_fake_vasp(wf, ir_ref_dirs, params_to_check=["ENCUT", "ISIF", "IBRION"])

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
        if mode not in ["H1-Ir_100 adsorbate optimization 1"]:
            raise ValueError("Invalid mode!")

        if "adsorbate" in mode:
            self.assertEqual(d["formula_reduced_abc"], "H1 Ir16")
        # Check relaxation of adsorbate
        # Check slab calculations
        # Check structure optimization

    def test_wf(self):
        self.wf_1 = self._simulate_vasprun(self.wf_1)

        self.assertEqual(len(self.wf_1.fws), 5)
        # check vasp parameters for ionic relaxation
        defo_vis = [fw.tasks[1]['vasp_input_set']
                    for fw in self.wf_1.fws if "adsorbate" in fw.name]
        assert all([vis.user_incar_settings['EDIFFG']==-0.05 for vis in defo_vis])
        assert all([vis.user_incar_settings['ISIF']==0 for vis in defo_vis])
        self.lp.add_wf(self.wf_1)
        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))

        # check relaxation
        d = self._get_task_collection().find_one({"task_label": "H1-Ir_100 adsorbate optimization 1"})
        self._check_run(d, mode="H1-Ir_100 adsorbate optimization 1")

        wf = self.lp.get_wf_by_fw_id(1)
        self.assertTrue(all([s == 'COMPLETED' for s in wf.fw_states.values()]))


if __name__ == "__main__":
    unittest.main()
