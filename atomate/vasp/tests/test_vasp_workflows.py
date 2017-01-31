# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import json
import os
import shutil
import unittest
import zlib

import gridfs
from pymongo import MongoClient, DESCENDING

from fireworks import LaunchPad, FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.vasp.powerups import use_custodian, add_namefile, use_fake_vasp, add_trackers, \
    add_bandgap_check
from atomate.vasp.workflows.base.core import get_wf

from pymatgen import SETTINGS
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.util.testing import PymatgenTest

__author__ = 'Anubhav Jain, Kiran Mathew'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "common", "reference_files", "db_connections")
reference_dir = os.path.join(module_dir, "reference_files")

ref_dirs_si = {"structure optimization": os.path.join(reference_dir, "Si_structure_optimization"),
             "static": os.path.join(reference_dir, "Si_static"),
             "nscf uniform": os.path.join(reference_dir, "Si_nscf_uniform"),
             "nscf line": os.path.join(reference_dir, "Si_nscf_line")}

DEBUG_MODE = True  # If true, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...


class TestVaspWorkflows(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # TODO: update this for the latest pymatgen...
        if not SETTINGS.get("PMG_VASP_PSP_DIR"):
            SETTINGS["PMG_VASP_PSP_DIR"] = os.path.join(module_dir, "reference_files")
            print('This system is not set up to run VASP jobs. '
                  'Please set PMG_VASP_PSP_DIR variable in your ~/.pmgrc.yaml file.')

        cls.struct_si = PymatgenTest.get_structure("Si")

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

    def _check_run(self, d, mode):
        if mode not in ["structure optimization", "static", "nscf uniform", "nscf line"]:
            raise ValueError("Invalid mode!")

        self.assertEqual(d["formula_pretty"], "Si")
        self.assertEqual(d["formula_anonymous"], "A")
        self.assertEqual(d["nelements"], 1)
        self.assertEqual(d["state"], "successful")
        self.assertAlmostEqual(d["calcs_reversed"][0]["output"]["structure"]["lattice"]["a"], 3.867, 2)
        self.assertEqual(d["output"]["is_gap_direct"], False)

        if mode in ["structure optimization", "static"]:
            self.assertAlmostEqual(d["output"]["energy"], -10.850, 2)
            self.assertAlmostEqual(d["output"]["energy_per_atom"], -5.425, 2)

        elif mode in ["ncsf uniform"]:
            self.assertAlmostEqual(d["output"]["energy"], -10.828, 2)
            self.assertAlmostEqual(d["output"]["energy_per_atom"], -5.414, 2)

        self.assertAlmostEqual(d["output"]["bandgap"], 0.65, 1)

        if "nscf" in mode:
            self.assertEqual(d["calcs_reversed"][0]["output"]["outcar"]["total_magnetization"], None)
        else:
            self.assertAlmostEqual(d["calcs_reversed"][0]["output"]["outcar"]["total_magnetization"], 0, 3)

        self.assertLess(d["run_stats"]["overall"]["Elapsed time (sec)"], 180)  # run should take under 3 minutes

        # check the DOS and band structure
        if mode == "nscf uniform" or mode == "nscf line":
            fs = gridfs.GridFS(self._get_task_database(), 'bandstructure_fs')

            # check the band structure
            bs_fs_id = d["calcs_reversed"][0]["bandstructure_fs_id"]
            bs_json = zlib.decompress(fs.get(bs_fs_id).read())
            bs = json.loads(bs_json.decode())
            self.assertEqual(bs["is_spin_polarized"], False)
            self.assertEqual(bs["band_gap"]["direct"], False)
            self.assertAlmostEqual(bs["band_gap"]["energy"], 0.65, 1)
            self.assertEqual(bs["is_metal"], False)

            if mode == "nscf uniform":
                for k in ["is_spin_polarized", "band_gap", "structure",
                          "kpoints", "is_metal", "vbm", "cbm", "labels_dict",
                          "projections", "lattice_rec", "bands"]:
                    self.assertTrue(k in bs)
                    self.assertIsNotNone(bs[k])

                self.assertEqual(bs["@class"], "BandStructure")

            else:
                for k in ["is_spin_polarized", "band_gap", "structure",
                          "kpoints", "is_metal", "vbm", "cbm", "labels_dict",
                          "projections", "lattice_rec", "bands", "branches"]:
                    self.assertTrue(k in bs)
                    self.assertIsNotNone(bs[k])
                self.assertEqual(bs["@class"], "BandStructureSymmLine")

            # check the DOS
            if mode == "nscf uniform":
                fs = gridfs.GridFS(self._get_task_database(), 'dos_fs')
                dos_fs_id = d["calcs_reversed"][0]["dos_fs_id"]

                dos_json = zlib.decompress(fs.get(dos_fs_id).read())
                dos = json.loads(dos_json.decode())
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
        structure = self.struct_si
        my_wf = get_wf(structure, "optimize_only.yaml", vis=MPRelaxSet(structure, force_gamma=True),
                       common_params={"vasp_cmd": VASP_CMD})
        if not VASP_CMD:
            my_wf = use_fake_vasp(my_wf, ref_dirs_si)
        else:
            my_wf = use_custodian(my_wf)
        self.lp.add_wf(my_wf)

        # run the workflow
        rapidfire(self.lp)

        fw = self.lp.get_fw_by_id(1)

        with open(os.path.join(fw.launches[-1].launch_dir, "task.json")) as f:
            d = json.load(f)
            self._check_run(d, mode="structure optimization")

    def test_single_Vasp_dbinsertion(self):
        # add the workflow
        structure = self.struct_si
        # instructs to use db_file set by FWorker, see env_chk
        my_wf = get_wf(structure, "optimize_only.yaml", vis=MPRelaxSet(structure, force_gamma=True),
                       common_params={"vasp_cmd": VASP_CMD,
                                      "db_file": ">>db_file<<"})
        if not VASP_CMD:
            my_wf = use_fake_vasp(my_wf, ref_dirs_si)
        else:
            my_wf = use_custodian(my_wf)
        self.lp.add_wf(my_wf)

        # run the workflow
        # set the db_file variable
        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))

        d = self._get_task_collection().find_one()
        self._check_run(d, mode="structure optimization")

    def test_bandstructure_Vasp(self):
        # add the workflow
        structure = self.struct_si
        # instructs to use db_file set by FWorker, see env_chk
        my_wf = get_wf(structure, "bandstructure.yaml",
                       vis=MPRelaxSet(structure, force_gamma=True),
                       common_params={"vasp_cmd": VASP_CMD,
                                      "db_file": ">>db_file<<"})
        if not VASP_CMD:
            my_wf = use_fake_vasp(my_wf, ref_dirs_si)
        else:
            my_wf = use_custodian(my_wf)

        my_wf = add_namefile(my_wf)  # add a slug of fw-name to output files

        self.lp.add_wf(my_wf)

        # run the workflow
        # set the db_file variable
        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))

        # make sure the structure relaxation ran OK
        d = self._get_task_collection().find_one({"task_label": "structure optimization"},
                                                 sort=[("_id", DESCENDING)])
        self._check_run(d, mode="structure optimization")

        # make sure the static run ran OK
        d = self._get_task_collection().find_one({"task_label": "static"}, sort=[("_id", DESCENDING)])
        self._check_run(d, mode="static")

        # make sure the uniform run ran OK
        d = self._get_task_collection().find_one({"task_label": "nscf uniform"}, sort=[("_id", DESCENDING)])
        self._check_run(d, mode="nscf uniform")

        # make sure the uniform run ran OK
        d = self._get_task_collection().find_one({"task_label": "nscf line"}, sort=[("_id", DESCENDING)])
        self._check_run(d, mode="nscf line")

    def test_bandgap_check_Vasp(self):
        # add the workflow
        structure = self.struct_si
        # instructs to use db_file set by FWorker, see env_chk
        my_wf = get_wf(structure, "bandstructure.yaml",
                       vis=MPRelaxSet(structure, force_gamma=True),
                       common_params={"vasp_cmd": VASP_CMD,
                                      "db_file": ">>db_file<<"})
        if not VASP_CMD:
            my_wf = use_fake_vasp(my_wf, ref_dirs_si)
        else:
            my_wf = use_custodian(my_wf)

        my_wf = add_namefile(my_wf)  # add a slug of fw-name to output files
        my_wf = add_bandgap_check(my_wf, check_bandgap_params={"max_gap": 0.1}, fw_name_constraint="structure optimization")
        self.lp.add_wf(my_wf)

        # run the workflow
        # set the db_file variable
        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))

        # structure optimization should be completed
        self.assertEqual(self.lp.fireworks.find_one(
            {"name": "Si-structure optimization"}, {"state": 1})["state"],
                         "COMPLETED")

        self.assertEqual(self.lp.fireworks.find_one(
            {"name": "Si-static"}, {"state": 1})["state"],
                         "DEFUSED")

    def test_trackers(self):
        # add the workflow
        structure = self.struct_si
        my_wf = get_wf(structure, "optimize_only.yaml",
                       vis=MPRelaxSet(structure, force_gamma=True),
                       common_params={"vasp_cmd": VASP_CMD})

        if not VASP_CMD:
            my_wf = use_fake_vasp(my_wf, ref_dirs_si)
        else:
            my_wf = use_custodian(my_wf)

        my_wf = add_trackers(my_wf)
        self.lp.add_wf(my_wf)

        # run the workflow
        rapidfire(self.lp)

        for x in self.lp.get_tracker_data(1):
            for t in x["trackers"]:
                self.assertGreater(len(t.content.split("\n")), 20)


if __name__ == "__main__":
    unittest.main()
