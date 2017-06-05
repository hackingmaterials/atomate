# coding: utf-8

from __future__ import unicode_literals

import unittest
from io import open
import os
import json
import shutil

from pymongo import MongoClient

from fireworks import LaunchPad

MODULE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)))
DB_DIR = os.path.join(MODULE_DIR, "..", "common", "test_files")

DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...


class AtomateTest(unittest.TestCase):

    def setUp(self):
        self.module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
        self.scratch_dir = os.path.join(self.module_dir, "scratch")
        if os.path.exists(self.scratch_dir):
            shutil.rmtree(self.scratch_dir)
        os.makedirs(self.scratch_dir)
        os.chdir(self.scratch_dir)
        try:
            self.lp = LaunchPad.from_file(os.path.join(DB_DIR, "my_launchpad.yaml"))
            self.lp.reset("", require_password=False)
        except:
            raise unittest.SkipTest('Cannot connect to MongoDB! Is the database server running? '
                                    'Are the credentials correct?')

    def get_task_database(self):
        with open(os.path.join(DB_DIR, "db.json")) as f:
            creds = json.loads(f.read())
            conn = MongoClient(creds["host"], creds["port"])
            db = conn[creds["database"]]
            if "admin_user" in creds:
                db.authenticate(creds["admin_user"], creds["admin_password"])
            return db

    def get_task_collection(self, coll_name=None):
        with open(os.path.join(DB_DIR, "db.json")) as f:
            creds = json.loads(f.read())
            db = self.get_task_database()
            coll_name = coll_name or creds["collection"]
            return db[coll_name]

    def tearDown(self):
        if not DEBUG_MODE:
            shutil.rmtree(self.scratch_dir)
            self.lp.reset("", require_password=False)
            db = self.get_task_database()
            for coll in db.collection_names():
                if coll != "system.indexes":
                    db[coll].drop()
            os.chdir(self.module_dir)