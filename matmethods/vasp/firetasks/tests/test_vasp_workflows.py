import json
import os
import shutil

from pymongo import MongoClient

from fireworks import LaunchPad, FWorker
from fireworks.core.rocket_launcher import rapidfire
from matmethods.vasp.examples.vasp_workflows import get_wf_single_Vasp, get_wf_double_Vasp
from matmethods.vasp.firetasks.tests.vasp_fake import make_fake_workflow
from pymatgen import IStructure, Lattice
from pymatgen.io.vasp.sets import MPVaspInputSet

import unittest

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "reference_files", "db_connections")
DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test

class TestVaspWorkflows(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if not os.environ.get("VASP_PSP_DIR"):
            raise unittest.SkipTest('This system is not set up to run VASP jobs. Please set your VASP_PSP_DIR environment variable.')
        try:
            cls.lp = LaunchPad.from_file(os.path.join(db_dir, "my_launchpad.yaml"))
        except:
            raise unittest.SkipTest('Cannot connect to MongoDB! Is the database server running? Are the credentials correct?')
        cls.lp.reset("", require_password=False)

        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                                [1.9200989668, 3.3257101909, 0.00],
                                [0.00, -2.2171384943, 3.1355090603]])
        cls.struct_si = IStructure(lattice, ["Si"] * 2, coords)

        cls.module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
        cls.scratch_dir = os.path.join(cls.module_dir, "scratch")

    def setUp(self):
        if os.path.exists(self.scratch_dir):
            shutil.rmtree(self.scratch_dir)
        os.makedirs(self.scratch_dir)
        os.chdir(self.scratch_dir)

    def tearDown(self):
        if not DEBUG_MODE:
            shutil.rmtree(self.scratch_dir)
            self._get_task_collection().delete_many({})
            self.lp.reset("", require_password=False)

    def _get_task_collection(self):
        with open(os.path.join(db_dir, "db.json")) as f:
            creds = json.loads(f.read())
            conn = MongoClient(creds["host"], creds["port"])
            db = conn[creds["database"]]
            if "readonly_user" in creds:
                db.authenticate(creds["readonly_user"], creds["readonly_password"])
            return db[creds["collection"]]

    def _check_relaxation_run(self, d):
        self.assertEqual(d["pretty_formula"], "Si")
        self.assertEqual(d["nelements"], 1)
        self.assertEqual(d["state"], "successful")
        self.assertAlmostEqual(d["output"]["final_energy"], -10.850, 2)
        self.assertAlmostEqual(d["output"]["final_energy_per_atom"], -5.425, 2)
        self.assertAlmostEqual(d["calculations"][0]["output"]["crystal"]["lattice"]["a"], 3.867, 2)
        self.assertAlmostEqual(d["calculations"][0]["output"]["outcar"]["total_magnetization"], 0, 3)
        self.assertAlmostEqual(d["analysis"]["bandgap"], 0.85, 1)
        self.assertEqual(d["analysis"]["is_gap_direct"], True)
        self.assertLess(d["run_stats"]["overall"]["Elapsed time (sec)"], 180)  # run should take under 3 minutes


    def test_single_Vasp(self):
        # add the workflow
        vis = MPVaspInputSet()
        structure = self.struct_si
        my_wf = get_wf_single_Vasp(structure, vis)
        my_wf = make_fake_workflow(my_wf)
        self.lp.add_wf(my_wf)

        # run the workflow
        rapidfire(self.lp)

        fw = self.lp.get_fw_by_id(1)

        with open(os.path.join(fw.launches[-1].launch_dir, "task.json")) as f:
            d = json.load(f)
            self._check_relaxation_run(d)

    def test_single_Vasp_dbinsertion(self):
        # add the workflow
        vis = MPVaspInputSet()
        structure = self.struct_si
        my_wf = get_wf_single_Vasp(structure, vis, db_file=">>db_file<<")  # instructs to use db_file set by FWorker, see env_chk
        my_wf = make_fake_workflow(my_wf)
        self.lp.add_wf(my_wf)

        # run the workflow
        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))  # set the db_file variable

        d = self._get_task_collection().find_one()
        self._check_relaxation_run(d)

    def test_double_Vasp(self):
        # add the workflow
        vis = MPVaspInputSet()
        structure = self.struct_si
        my_wf = get_wf_double_Vasp(structure, vis)
        my_wf = make_fake_workflow(my_wf)
        self.lp.add_wf(my_wf)

        # run the workflow
        rapidfire(self.lp)

        # make sure the structure relaxation ran OK
        fw_id = self.lp.get_fw_ids({"name": "structure optimization"})[0]
        fw = self.lp.get_fw_by_id(fw_id)
        with open(os.path.join(fw.launches[-1].launch_dir, "task.json")) as f:
            d = json.load(f)
            self._check_relaxation_run(d)
