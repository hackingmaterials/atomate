# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import json
import os
import shutil
import unittest

import numpy as np

from pymongo import MongoClient

from pymatgen import Structure
from pymatgen.io.feff.inputs import Tags

from fireworks.core.fworker import FWorker
from fireworks.core.launchpad import LaunchPad
from fireworks.core.rocket_launcher import rapidfire

from atomate.feff.workflows.xas import get_wf_xas


__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "common", "reference_files", "db_connections")
DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test
FEFF_CMD = None  # "feff"


class TestXASWorkflow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # CoO
        cls.structure = Structure.from_file(os.path.join(module_dir, "reference_files", "Co2O2.cif"))
        #PymatgenTest.get_mp_structure("mp-715460")
        cls.user_tag_settings = {"RPATH": -1,
                                 "SCF": "7 0 30 0.2 3",
                                 "FMS": "9 0",
                                 "LDOS": "-30.0 30.0 0.1",
                                 "RECIPROCAL": "",
                                 "EDGE": "K"}
        # 3rd site
        cls.absorbing_atom = 2
        cls.edge = "K"
        cls.nkpts = 1000
        cls.scratch_dir = os.path.join(module_dir, "scratch")

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

    def test_xas_wflow_abatom_by_idx(self):
        if not FEFF_CMD:
            # fake run
            feff_bin = "cp  ../../reference_files/xmu.dat ."
        else:
            feff_bin = FEFF_CMD

        wf = get_wf_xas(self.absorbing_atom, self.structure, spectrum_type="XANES", edge="K",
                        feff_cmd=feff_bin, db_file=">>db_file<<", use_primitive=False,
                        user_tag_settings=self.user_tag_settings)
        self.assertEqual(len(wf.as_dict()["fws"]), 1)

        self.lp.add_wf(wf)
        # run
        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))

        d = self._get_task_collection().find_one({"spectrum_type": "XANES"})
        self._check_run(d)

    def test_xas_wflow_abatom_by_symbol(self):
        wf_prim = get_wf_xas("O", self.structure, spectrum_type="XANES", edge="K",
                             use_primitive=True, user_tag_settings=self.user_tag_settings)
        wf = get_wf_xas("O", self.structure, spectrum_type="XANES", edge="K",
                        use_primitive=False, user_tag_settings=self.user_tag_settings)
        self.assertEqual(len(wf_prim.as_dict()["fws"]), 1)
        self.assertEqual(len(wf.as_dict()["fws"]), 2)

    def test_xanes_vs_exafs(self):
        wf_xanes = get_wf_xas(self.absorbing_atom, self.structure, spectrum_type="XANES", edge="K",
                              user_tag_settings=self.user_tag_settings)
        wf_exafs = get_wf_xas(self.absorbing_atom, self.structure, spectrum_type="EXAFS", edge="K")

        self.assertEqual(wf_xanes.as_dict()["fws"][0]["spec"]['_tasks'][0]['feff_input_set']['@class'],
                         'MPXANESSet')
        self.assertEqual(wf_exafs.as_dict()["fws"][0]["spec"]['_tasks'][0]['feff_input_set']['@class'],
                         'MPEXAFSSet')

    def _check_run(self, d):
        run_dir = d["dir_name"]
        self.assertEqual(d["edge"], self.edge)
        self.assertEqual(d["absorbing_atom"], self.absorbing_atom)
        tags = Tags.from_file(os.path.join(run_dir, "feff.inp"))
        self.assertEqual(d["input_parameters"], tags.as_dict())
        xmu = np.loadtxt((os.path.join(run_dir, "xmu.dat")))
        self.assertEqual(d["spectrum"], xmu.tolist())

    def _get_task_database(self):
        with open(os.path.join(db_dir, "db.json")) as f:
            creds = json.loads(f.read())
            conn = MongoClient(creds["host"], creds["port"])
            db = conn[creds["database"]]
            if "admin_user" in creds:
                db.authenticate(creds["admin_user"], creds["admin_password"])
            return db

    def _get_task_collection(self):
        with open(os.path.join(db_dir, "db.json")) as f:
            creds = json.loads(f.read())
            db = self._get_task_database()
            return db[creds["collection"]]

    def tearDown(self):
        if not DEBUG_MODE:
            shutil.rmtree(self.scratch_dir)
            self.lp.reset("", require_password=False)
            db = self._get_task_database()
            for coll in db.collection_names():
                if coll != "system.indexes":
                    db[coll].drop()


if __name__ == "__main__":
    unittest.main()
