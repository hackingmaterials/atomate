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

from atomate.utils.testing import AtomateTest
from atomate.vasp.powerups import use_fake_vasp
from atomate.vasp.firetasks.parse_outputs import PolarizationToDb

from pymatgen import SETTINGS

__author__ = 'Tess Smidt'
__email__ = 'blondegeek@gmail.com'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
ref_dir = os.path.join(module_dir, "..", "..", "test_files")

from pymatgen.core.structure import Structure

DEBUG_MODE = True  # If true, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...


class TestFerroelectricWorkflow(AtomateTest):

    def test_polarizationtodb(self):

        import bson
        import gzip

        reference_dir = os.path.abspath(os.path.join(ref_dir, "ferroelectric_wf"))

        with gzip.open(os.path.join(reference_dir, "tasks.bson.gz")) as f:
            coll_raw = f.read()

        coll = bson.decode_all(coll_raw)

        db = self.get_task_collection()
        for c in coll:
            db.insert(c)

        new_fw_spec = {'_fw_env': {"db_file": os.path.join(db_dir, "db.json")},
                       'tags':['wfid_1494203093.06934658']}

        analysis = PolarizationToDb(db_file='>>db_file<<', name="_polarization_post_processing")
        analysis.run_task(new_fw_spec)

        # Check recovered change in polarization
        coll = self.get_task_collection("polarization_tasks")
        d = coll.find_one()
        self.assertAlmostEqual(d['polarization_change_norm'], 46.288752795325244, 5)



if __name__ == "__main__":
    unittest.main()
