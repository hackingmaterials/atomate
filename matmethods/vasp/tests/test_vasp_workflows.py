import json
import os
import shutil

import gridfs
import zlib
from pymongo import MongoClient

from fireworks import LaunchPad, FWorker
from fireworks.core.rocket_launcher import rapidfire
from matmethods.vasp.examples.vasp_workflows import get_wf_single_Vasp, get_wf_bandstructure_Vasp
from matmethods.vasp.firetasks.new_input_sets import StructureOptimizationVaspInputSet
from matmethods.vasp.tests.vasp_fake import make_fake_workflow, make_custodian_workflow
from pymatgen import IStructure, Lattice

import unittest

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "reference_files", "db_connections")
DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...

class TestVaspWorkflows(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if not os.environ.get("VASP_PSP_DIR"):
            raise unittest.SkipTest('This system is not set up to run VASP jobs. Please set your VASP_PSP_DIR environment variable.')

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
        try:
            self.lp = LaunchPad.from_file(os.path.join(db_dir, "my_launchpad.yaml"))
        except:
            raise unittest.SkipTest('Cannot connect to MongoDB! Is the database server running? Are the credentials correct?')
        self.lp.reset("", require_password=False)

    def tearDown(self):
        if not DEBUG_MODE:
            shutil.rmtree(self.scratch_dir)
            self.lp.reset("", require_password=False)
            db = self._get_task_database()
            for coll in db.collection_names():
                if coll != "system.indexes":
                    db[coll].drop()

    def _get_task_database(self):
        with open(os.path.join(db_dir, "db.json")) as f:
            creds = json.loads(f.read())
            conn = MongoClient(creds["host"], creds["port"])
            db = conn[creds["database"]]
            if "readonly_user" in creds:
                db.authenticate(creds["readonly_user"], creds["readonly_password"])
            return db

    def _get_task_collection(self):
        with open(os.path.join(db_dir, "db.json")) as f:
            creds = json.loads(f.read())
            db = self._get_task_database()
            return db[creds["collection"]]

    def _check_run(self, d, mode):
        if mode not in ["structure optimization", "static", "nscf uniform", "nscf line"]:
            raise ValueError("Invalid mode!")

        self.assertEqual(d["pretty_formula"], "Si")
        self.assertEqual(d["nelements"], 1)
        self.assertEqual(d["state"], "successful")
        self.assertAlmostEqual(d["calculations"][0]["output"]["crystal"]["lattice"]["a"], 3.867, 2)
        self.assertEqual(d["analysis"]["is_gap_direct"], False)

        if mode != "nscf line":
            self.assertAlmostEqual(d["output"]["final_energy"], -10.850, 2)
            self.assertAlmostEqual(d["output"]["final_energy_per_atom"], -5.425, 2)

        if "nscf" in mode:
            self.assertEqual(d["calculations"][0]["output"]["outcar"]["total_magnetization"], None)
            self.assertAlmostEqual(d["analysis"]["bandgap"], 0.62, 1)
        else:
            self.assertAlmostEqual(d["calculations"][0]["output"]["outcar"]["total_magnetization"], 0, 3)
            self.assertAlmostEqual(d["analysis"]["bandgap"], 0.65, 1)

        self.assertLess(d["run_stats"]["overall"]["Elapsed time (sec)"], 180)  # run should take under 3 minutes

        # check the DOS and band structure
        if mode == "nscf uniform" or mode == "nscf line":
            fs = gridfs.GridFS(self._get_task_database(), 'bandstructure_fs')

            # check the band structure
            bs_fs_id = d["calculations"][0]["bandstructure_fs_id"]
            bs_json = zlib.decompress(fs.get(bs_fs_id).read())
            bs = json.loads(bs_json)

            self.assertEqual(bs["is_spin_polarized"], True)  # TODO: should this be false?
            self.assertEqual(bs["band_gap"]["direct"], True)  # TODO: this should almost certainly be false
            self.assertAlmostEqual(bs["band_gap"]["energy"], 0.85, 1)  # TODO: this does not match the value from earlier in the unit test for nscf uniform runs
            self.assertEqual(bs["is_metal"], False)

            if mode == "nscf uniform":
                for k in ["is_spin_polarized", "band_gap", "structure", "kpoints", "is_metal", "vbm", "cbm", "labels_dict", "projections", "lattice_rec", "bands"]:
                    self.assertTrue(k in bs)
                    self.assertIsNotNone(bs[k])

                self.assertEqual(bs["@class"], "BandStructure")

            else:
                for k in ["is_spin_polarized", "band_gap", "structure", "kpoints", "is_metal", "vbm", "cbm", "labels_dict", "projections", "lattice_rec", "bands", "branches"]:
                    self.assertTrue(k in bs)
                    self.assertIsNotNone(bs[k])
                self.assertEqual(bs["@class"], "BandStructureSymmLine")
                # TODO: the branches key seems wrong??

            # check the DOS
            if mode == "nscf uniform":
                fs = gridfs.GridFS(self._get_task_database(), 'dos_fs')
                dos_fs_id = d["calculations"][0]["dos_fs_id"]

                dos_json = zlib.decompress(fs.get(dos_fs_id).read())
                dos = json.loads(dos_json)
                for k in ["densities", "energies", "pdos", "spd_dos", "atom_dos", "structure"]:
                    self.assertTrue(k in dos)
                    self.assertIsNotNone(dos[k])

                self.assertAlmostEqual(dos["spd_dos"]["p"]["efermi"], 5.625, 1)
                self.assertAlmostEqual(dos["atom_dos"]["Si"]["efermi"], 5.625, 1)
                self.assertAlmostEqual(dos["structure"]["lattice"]["a"], 3.867, 2)
                self.assertAlmostEqual(dos["spd_dos"]["p"]["efermi"], 5.625, 1)
                self.assertAlmostEqual(dos["atom_dos"]["Si"]["efermi"], 5.625, 1)
                self.assertAlmostEqual(dos["structure"]["lattice"]["a"], 3.867, 2)

    def test_single_Vasp(self):
        # add the workflow
        vis = StructureOptimizationVaspInputSet()
        structure = self.struct_si
        my_wf = get_wf_single_Vasp(structure, vis, vasp_cmd=VASP_CMD, task_label="structure optimization")
        if not VASP_CMD:
            my_wf = make_fake_workflow(my_wf)
        else:
            my_wf = make_custodian_workflow(my_wf)
        self.lp.add_wf(my_wf)

        # run the workflow
        rapidfire(self.lp)

        fw = self.lp.get_fw_by_id(1)

        with open(os.path.join(fw.launches[-1].launch_dir, "task.json")) as f:
            d = json.load(f)
            self._check_run(d, mode="structure optimization")

    def test_single_Vasp_dbinsertion(self):
        # add the workflow
        vis = StructureOptimizationVaspInputSet()
        structure = self.struct_si
        my_wf = get_wf_single_Vasp(structure, vis, vasp_cmd=VASP_CMD, db_file=">>db_file<<", task_label="structure optimization")  # instructs to use db_file set by FWorker, see env_chk
        if not VASP_CMD:
            my_wf = make_fake_workflow(my_wf)
        else:
            my_wf = make_custodian_workflow(my_wf)
        self.lp.add_wf(my_wf)

        # run the workflow
        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))  # set the db_file variable

        d = self._get_task_collection().find_one()
        self._check_run(d, mode="structure optimization")

    def test_bandstructure_Vasp(self):
        # add the workflow
        vis = StructureOptimizationVaspInputSet()
        structure = self.struct_si
        my_wf = get_wf_bandstructure_Vasp(structure, vis, vasp_cmd=VASP_CMD, db_file=">>db_file<<")  # instructs to use db_file set by FWorker, see env_chk
        if not VASP_CMD:
            my_wf = make_fake_workflow(my_wf)
        else:
            my_wf = make_custodian_workflow(my_wf)
        self.lp.add_wf(my_wf)

        # run the workflow
        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))  # set the db_file variable

        # make sure the structure relaxation ran OK
        d = self._get_task_collection().find_one({"task_label": "structure optimization"})
        self._check_run(d, mode="structure optimization")

        # make sure the static run ran OK
        d = self._get_task_collection().find_one({"task_label": "static"})
        self._check_run(d, mode="static")

        # make sure the uniform run ran OK
        d = self._get_task_collection().find_one({"task_label": "nscf uniform"})
        self._check_run(d, mode="nscf uniform")

        # make sure the uniform run ran OK
        d = self._get_task_collection().find_one({"task_label": "nscf line"})
        self._check_run(d, mode="nscf line")

if __name__ == "__main__":
    unittest.main()