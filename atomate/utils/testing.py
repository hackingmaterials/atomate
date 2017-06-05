# coding: utf-8

from __future__ import unicode_literals

import unittest
from io import open
import os
import json

from pymongo import MongoClient

from pymatgen import SETTINGS

if not SETTINGS.get("PMG_VASP_PSP_DIR"):
    SETTINGS["PMG_VASP_PSP_DIR"] = os.path.join(module_dir, "..", "..", "tests", "reference_files")
    print('This system is not set up to run VASP jobs. '
          'Please set PMG_VASP_PSP_DIR variable in your ~/.pmgrc.yaml file.')

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "common", "test_files")


class AtomateTest(unittest.TestCase):

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
