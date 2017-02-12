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
from atomate.vasp.workflows.presets.core import wf_gibbs_free_energy

from pymatgen import SETTINGS
from pymatgen.util.testing import PymatgenTest
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

__author__ = 'Brandon Bocklund'
__email__ = 'brandonbocklund@gmail.com'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "reference_files", "db_connections")
ref_dir = os.path.join(module_dir, "test_files")

DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...

class TestGibbsWorkflowPhonon(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not SETTINGS.get("PMG_VASP_PSP_DIR"):
            SETTINGS["PMG_VASP_PSP_DIR"] = os.path.join(module_dir, "..", "..", "tests", "reference_files")
            print('This system is not set up to run VASP jobs. '
                  'Please set PMG_VASP_PSP_DIR variable in your ~/.pmgrc.yaml file.')

        cls.struct_si = SpacegroupAnalyzer(
                PymatgenTest.get_structure("Si")).get_conventional_standard_structure()
        cls.struct_si.make_supercell(2)
        cls.scratch_dir = os.path.join(module_dir, "scratch")
        cls.gibbs_config = {"deformations":[(np.identity(3)*(1+x)**(1.0/3.0)).tolist() for x in np.linspace(-0.05, 0.05, 5)],
                            "t_min":0, "t_step": 10, "t_max":1500,
                            "mesh": [5, 5, 5], # small mesh for faster tests
                            "qha_type":"phonopy",
                            "vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<"}
        cls.wf = wf_gibbs_free_energy(cls.struct_si, c=cls.gibbs_config)

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
        reference_dir = os.path.abspath(os.path.join(ref_dir, "gibbs_wf_phonon"))
        si_ref_dirs = {"structure optimization": os.path.join(reference_dir, "1"),
                       "gibbs deformation 0": os.path.join(reference_dir, "2"),
                       "gibbs deformation 1": os.path.join(reference_dir, "3"),
                       "gibbs deformation 2": os.path.join(reference_dir, "4"),
                       "gibbs deformation 3": os.path.join(reference_dir, "5"),
                       "gibbs deformation 4": os.path.join(reference_dir, "6")}
        return use_fake_vasp(wf, si_ref_dirs, params_to_check=["LREAL", "ADDGRID", "IBRION"])

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
        if mode not in ["structure optimization", "gibbs deformation 0",
                        "gibbs deformation 2", "gibbs analysis"]:
            raise ValueError("Invalid mode!")

        self.assertEqual(d["formula_pretty"], "Si")
        self.assertEqual(d["nsites"], 64)

        if mode not in ["gibbs analysis"]:
            self.assertEqual(d["formula_anonymous"], "A")
            self.assertEqual(d["nelements"], 1)

        if mode in ["structure optimization"]:
            self.assertAlmostEqual(d["calcs_reversed"][-1]["output"]["structure"]["lattice"]["a"], 10.937, 2)
            self.assertAlmostEqual(d["output"]["energy_per_atom"], -5.425, 2)

        elif mode in ["gibbs deformation 0"]:
            forces_5_0 = np.array([[-5.711e-05, 0.0, 0.0], [0.0, -0.01412188, -0.0], [0.0, -0.0, -5.711e-05]]) # check that the output is correct
            np.testing.assert_allclose(d["calcs_reversed"][0]["output"]['force_constants'][5][0], forces_5_0, rtol=1e-2)
            self.assertAlmostEqual(d["output"]["energy_per_atom"], -5.409, 2)
            self.assertAlmostEqual(d["calcs_reversed"][-1]["output"]["structure"]["lattice"]["a"], 10.752, 2)

        elif mode in ["gibbs deformation 2"]:
            forces_last_0 = np.array([[-0.00056991, 0.00086526, -0.00086526], [0.00086526, -0.00056991, -0.00086526], [-0.00086526, -0.00086526, -0.00056991]])
            np.testing.assert_allclose(d["calcs_reversed"][0]["output"]['force_constants'][-1][0], forces_last_0, rtol=1e-2)
            self.assertAlmostEqual(d["calcs_reversed"][-1]["output"]["structure"]["lattice"]["a"], 10.937, 2)
            self.assertAlmostEqual(d["output"]["energy_per_atom"], -5.425, 2)

        elif mode in ["gibbs analysis"]:
            self.assertEqual(len(d["force_constants"]),5)
            self.assertEqual(len(d["force_constants"]),len(d["energies"]))
            self.assertEqual(len(d["energies"]), len(d["volumes"]))
            self.assertEqual(len(d["gibbs_free_energy"]), 151)
            self.assertEqual(len(d["gibbs_free_energy"]), len(d["temperatures"]))
            self.assertEqual(d["temperatures"][-1], 1500, 'Max temperature differs from configured temperature')
            self.assertEqual(d["qha_type"], "phonopy")
            # we know that force constants, energies, and volumes are being calculated correctly so
            # testing of the calculated values should be handled by unit tests for analysis

    def test_wf(self):
        self.wf = self._simulate_vasprun(self.wf)
        self.assertEqual(len(self.wf.fws), 7)

        self.lp.add_wf(self.wf)

        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))

        reg = lambda name: {'$regex': name}
        # check relaxation
        d = self._get_task_collection().find_one({"task_label": reg("structure optimization")})
        self._check_run(d, mode="structure optimization")
        # check two of the deformation calculations
        d = self._get_task_collection().find_one({"task_label": reg("gibbs deformation 0")})
        self._check_run(d, mode="gibbs deformation 0")

        d = self._get_task_collection().find_one({"task_label": reg("gibbs deformation 2")})
        self._check_run(d, mode="gibbs deformation 2")

        # check the final results
        d = self._get_task_collection(coll_name="gibbs").find_one()
        if d is None: raise ValueError("Gibbs d not found")
        self._check_run(d, mode="gibbs analysis")


if __name__ == "__main__":
    unittest.main()
