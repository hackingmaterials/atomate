# coding: utf-8

import json
import os
import unittest
import shutil
from monty.os.path import which
from atomate.qchem.firetasks.critic2 import RunCritic2, ProcessCritic2
from atomate.utils.testing import AtomateTest
from custodian.qchem.handlers import QChemErrorHandler
from custodian.qchem.jobs import QCJob
from pymatgen.io.qchem.outputs import QCOutput
import numpy as np

__author__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

@unittest.skipIf(not which("critic2"), "critic2 executable not present")
class TestRunCritic2(AtomateTest):
    def setUp(self, lpad=False):
        os.chdir(os.path.join(module_dir, "..", "..", "test_files",
                                "critic_test_files", "small_critic_example"))
        out_file = "mol.qout"
        qc_out = QCOutput(filename=out_file)
        self.mol = qc_out.data["initial_molecule"]
        self.cube_file = "dens.0.cube"
        super(TestRunCritic2, self).setUp(lpad=False)

    def tearDown(self):
        os.remove("cpreport.json")
        os.remove("yt.json")

    def test_RunCritic2(self):
        os.chdir(os.path.join(module_dir, "..", "..", "test_files",
                                "critic_test_files", "small_critic_example"))
        firetask = RunCritic2(
            molecule=self.mol,
            cube_file="dens.0.cube.gz")
        firetask.run_task(fw_spec={})
        with open("cpreport_correct.json") as f:
            cpreport_reference = json.load(f)
        with open("yt_correct.json") as f:
            yt_reference = json.load(f)
        with open("cpreport.json") as f:
            cpreport = json.load(f)
        with open("yt.json") as f:
            yt = json.load(f)
        # Context for below - reference files were built before units were added
        # to Critic2, and we avoid testing the actual critical points because they
        # can change order between runs. But comparing everything else is sufficient.
        for key in cpreport:
            if key in ["structure", "field"]:
                self.assertEqual(cpreport_reference[key],cpreport[key])
        for key in yt:
            if key != "units":
                self.assertEqual(yt_reference[key],yt[key])


class TestProcessCritic2(AtomateTest):
    def setUp(self, lpad=False):
        os.chdir(os.path.join(module_dir, "..", "..", "test_files",
                                "critic_test_files", "critic_example"))
        out_file = "mol.qout"
        qc_out = QCOutput(filename=out_file)
        self.mol = qc_out.data["initial_molecule"]
        self.cube_file = "dens.0.cube.gz"
        shutil.move("bonding.json","bonding_correct.json")
        shutil.move("processed_critic2.json","processed_correct.json")
        super(TestProcessCritic2, self).setUp(lpad=False)

    def tearDown(self):
        os.remove("bonding.json")
        shutil.move("bonding_correct.json","bonding.json")
        os.remove("processed_critic2.json")
        shutil.move("processed_correct.json","processed_critic2.json")

    def test_ProcessCritic2(self):
        os.chdir(os.path.join(module_dir, "..", "..", "test_files",
                                "critic_test_files", "critic_example"))
        firetask = ProcessCritic2(
            molecule=self.mol)
        firetask.run_task(fw_spec={})
        with open("bonding_correct.json") as f:
            reference = json.load(f)
        with open("bonding.json") as f:
            just_built = json.load(f)
        self.assertEqual(reference,just_built)
        with open("processed_correct.json") as f:
            reference = json.load(f)
        with open("processed_critic2.json") as f:
            just_built = json.load(f)
        self.assertEqual(reference,just_built)


if __name__ == "__main__":
    unittest.main()
