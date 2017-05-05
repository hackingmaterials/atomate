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
from atomate.vasp.workflows.base.ferroelectric import get_wf_ferroelectric, get_wf_id

from pymatgen import SETTINGS

__author__ = 'Tess Smidt'
__email__ = 'blondegeek@gmail.com'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
ref_dir = os.path.join(module_dir, "test_files")

from pymatgen.core.structure import Structure

DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...

class TestFerroelectricWorkflow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not SETTINGS.get("PMG_VASP_PSP_DIR"):
            SETTINGS["PMG_VASP_PSP_DIR"] = os.path.join(module_dir, "..", "..", "tests", "reference_files")
            print('This system is not set up to run VASP jobs. '
                  'Please set PMG_VASP_PSP_DIR variable in your ~/.pmgrc.yaml file.')

        cls.bto_polar = Structure.from_file(ref_dir+"/ferroelectric_wf/"+"BTO_polar_POSCAR")
        cls.bto_nonpolar = Structure.from_file(ref_dir+"/ferroelectric_wf/"+"BTO_nonpolar_POSCAR")

        cls.wfid = "wfid_" + get_wf_id()

        cls.ferroelectric_config = {'vasp_cmd': '>>vasp_cmd<<',
                                    'db_file': '>>db_file<<',
                                    'nimages': 9,
                                    'relax' : True,
                                    'wfid': cls.wfid,
                                    'add_analysis_task': True}
        cls.wf = get_wf_ferroelectric(cls.bto_polar, cls.bto_nonpolar, **cls.ferroelectric_config)

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
        pass
        reference_dir = os.path.abspath(os.path.join(ref_dir, "ferroelectric_wf"))
        bto_ref_dirs = {"polar_relaxation": os.path.join(reference_dir, "polar_relaxation"),
                        "polar_static": os.path.join(reference_dir, "polar_static"),
                        "polar_polarization": os.path.join(reference_dir, "polar_polarization"),
                        "nonpolar_relaxation": os.path.join(reference_dir, "nonpolar_relaxation"),
                        "nonpolar_static": os.path.join(reference_dir, "nonpolar_static"),
                        "nonpolar_polarization": os.path.join(reference_dir, "nonpolar_polarization"),
                        "interpolation_1_static": os.path.join(reference_dir, "interpolation_1_static"),
                        "interpolation_2_static": os.path.join(reference_dir, "interpolation_2_static"),
                        "interpolation_3_static": os.path.join(reference_dir, "interpolation_3_static"),
                        "interpolation_4_static": os.path.join(reference_dir, "interpolation_4_static"),
                        "interpolation_5_static": os.path.join(reference_dir, "interpolation_5_static"),
                        "interpolation_6_static": os.path.join(reference_dir, "interpolation_6_static"),
                        "interpolation_7_static": os.path.join(reference_dir, "interpolation_7_static"),
                        "interpolation_8_static": os.path.join(reference_dir, "interpolation_8_static"),
                        "interpolation_1_polarization": os.path.join(reference_dir, "interpolation_1_polarization"),
                        "interpolation_2_polarization": os.path.join(reference_dir, "interpolation_2_polarization"),
                        "interpolation_3_polarization": os.path.join(reference_dir, "interpolation_3_polarization"),
                        "interpolation_4_polarization": os.path.join(reference_dir, "interpolation_4_polarization"),
                        "interpolation_5_polarization": os.path.join(reference_dir, "interpolation_5_polarization"),
                        "interpolation_6_polarization": os.path.join(reference_dir, "interpolation_6_polarization"),
                        "interpolation_7_polarization": os.path.join(reference_dir, "interpolation_7_polarization"),
                        "interpolation_8_polarization": os.path.join(reference_dir, "interpolation_8_polarization"),
                        "polarization_analysis": os.path.join(reference_dir, "polarization_analysis")}
        # Add test with for analysis?
        # Which params_to_check?
        return use_fake_vasp(wf, bto_ref_dirs, params_to_check=["ENCUT"])

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

        # Check polar and nonpolar relaxations
        if mode is 'polar_relaxation':
            self.assertAlmostEqual(d["calcs_reversed"][0]["output"]["structure"]["lattice"]["c"], 4.2157, 2)

        if mode is 'nonpolar_relaxation':
            self.assertAlmostEqual(d["calcs_reversed"][0]["output"]["structure"]["lattice"]["c"], 4.0330, 2)

        # Check interpolated structures
        if mode is 'interpolation_1_polarization':
            self.assertAlmostEqual(d["calcs_reversed"][0]["output"]["structure"]["lattice"]["c"], 4.1954, 2)

        if mode is 'interpolation_5_polarization':
            self.assertAlmostEqual(d["calcs_reversed"][0]["output"]["structure"]["lattice"]["c"], 4.1142, 2)

        # Check that Outcar has needed keys for polarization analysis.
        if '_polarization' in mode:

            # Check that Outcar has p_ion, p_elec, zval_dict
            assert d["calcs_reversed"][0]["output"]["outcar"].get("p_ion", None) is not None
            assert d["calcs_reversed"][0]["output"]["outcar"].get("p_elec", None) is not None
            assert d["calcs_reversed"][0]["output"]["outcar"].get("zval_dict", None) is not None

        # Check analysis step.
        if mode is "polarization_post_processing":
            self.assertAlmostEqual(d['polarization_change_norm'], 46.288752795325244)

    def test_wf(self):

        self.wf = self._simulate_vasprun(self.wf)

        # 2*relax + 10*polarization +1*polarization_analysis = 13
        self.assertEqual(len(self.wf.fws), 13)

        # check VASP parameters on polarization calculation for interpolated structures
        interpolated_polarization_vis = [fw.spec["_tasks"][8]['other_params']['lcalcpol']
            for fw in self.wf.fws if "polarization" in fw.name and "interpolation" in fw.name]

        assert all(interpolated_polarization_vis)

        self.lp.add_wf(self.wf)
        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))

        # Check polar relaxation
        d = self._get_task_collection().find_one({"task_label": "polar_relaxation"})
        self._check_run(d, "polar_relaxation")

        # Check nonpolar relaxation
        d = self._get_task_collection().find_one({"task_label": "nonpolar_relaxation"})
        self._check_run(d, "nonpolar_relaxation")

        # Check polarization calculations
        D = self._get_task_collection().find({"task_label": {"$regex": ".*polarization"}})
        for d in D:
            self._check_run(d, d["task_label"])

        # Check recovered change in polarization
        coll = self._get_task_collection("polarization_tasks")
        d = coll.find_one()
        self._check_run(d,"polarization_post_processing")

if __name__ == "__main__":
    unittest.main()
