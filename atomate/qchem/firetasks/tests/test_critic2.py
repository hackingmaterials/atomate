# coding: utf-8

import json
import os
import unittest
import shutil

from atomate.qchem.firetasks.critic2 import ProcessCritic2
from atomate.utils.testing import AtomateTest
from custodian.qchem.handlers import QChemErrorHandler
from custodian.qchem.jobs import QCJob
from pymatgen.io.qchem.outputs import QCOutput
import numpy as np

__author__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class TestProcessCritic2(AtomateTest):
    def setUp(self, lpad=False):
        os.chdir(os.path.join(module_dir, "..", "..", "test_files",
                                "critic_test_files", "critic_example"))
        out_file = "mol.qout"
        qc_out = QCOutput(filename=out_file)
        self.mol = qc_out.data["initial_molecule"]
        self.cube_file = "dens.0.cube.gz"
        shutil.move("bonding.json","bonding_correct.json")
        super(TestProcessCritic2, self).setUp(lpad=False)

    def tearDown(self):
        os.remove("bonding.json")
        shutil.move("bonding_correct.json","bonding.json")

    def test_ProcessCritic2(self):
        os.chdir(os.path.join(module_dir, "..", "..", "test_files",
                                "critic_test_files", "critic_example"))
        firetask = ProcessCritic2(
            molecule=self.mol,
            cp_name="CP.json",
            yt_name="YT.json")
        print(os.getcwd())
        firetask.run_task(fw_spec={})
        with open("bonding_correct.json") as f:
            reference = json.load(f)
        with open("bonding.json") as f:
            just_built = json.load(f)
        self.assertEqual(reference,just_built)


if __name__ == "__main__":
    unittest.main()
