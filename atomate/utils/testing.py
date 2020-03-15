# coding: utf-8


import unittest
from io import open
import os
import json
import shutil

from pymongo import MongoClient

from fireworks import LaunchPad

from pymatgen import SETTINGS

__author__ = 'Kiran Mathew'
__credits__ = 'Anubhav Jain'
__email__ = 'kmathew@lbl.gov'

MODULE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)))
DB_DIR = os.path.join(MODULE_DIR, "..", "common", "test_files")

DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...


class AtomateTest(unittest.TestCase):

    def setUp(self, lpad=True):
        """
        Create scratch directory(removes the old one if there is one) and change to it.
        Also initialize launchpad.
        """
        if not SETTINGS.get("PMG_VASP_PSP_DIR"):
            SETTINGS["PMG_VASP_PSP_DIR"] = os.path.abspath(os.path.join(MODULE_DIR, "..", "vasp", "test_files"))
            print('This system is not set up to run VASP jobs. '
                  'Please set PMG_VASP_PSP_DIR variable in your ~/.pmgrc.yaml file.')

        self.scratch_dir = os.path.join(MODULE_DIR, "scratch")
        if os.path.exists(self.scratch_dir):
            shutil.rmtree(self.scratch_dir)
        os.makedirs(self.scratch_dir)
        os.chdir(self.scratch_dir)
        if lpad:
            try:
                self.lp = LaunchPad.from_file(os.path.join(DB_DIR, "my_launchpad.yaml"))
                self.lp.reset("", require_password=False)
            except:
                raise unittest.SkipTest('Cannot connect to MongoDB! Is the database server running? '
                                        'Are the credentials correct?')

    # Note: the functions in matgendb.util, get_database and get_collection require db authentication
    # but the db.json config file used for atomate testing purpose doesnt require db authentication.
    # Hence the following 2 methods.
    def get_task_database(self):
        """
        Returns pymongo db connection.
        """
        with open(os.path.join(DB_DIR, "db.json")) as f:
            creds = json.loads(f.read())
            conn = MongoClient(creds["host"], creds["port"])
            db = conn[creds["database"]]
            if "admin_user" in creds:
                db.authenticate(creds["admin_user"], creds["admin_password"])
            return db

    def get_task_collection(self, coll_name=None):
        """
        Returns pymongo collection
        """
        with open(os.path.join(DB_DIR, "db.json")) as f:
            creds = json.loads(f.read())
            db = self.get_task_database()
            coll_name = coll_name or creds["collection"]
            return db[coll_name]

    def tearDown(self):
        """
        Remove the scratch directory and teardown the test db.
        """
        if not DEBUG_MODE:
            shutil.rmtree(self.scratch_dir)
            self.lp.reset("", require_password=False)
            db = self.get_task_database()
            for coll in db.collection_names():
                if coll != "system.indexes":
                    db[coll].drop()
            os.chdir(MODULE_DIR)
