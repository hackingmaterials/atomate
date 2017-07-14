# coding: utf-8
from __future__ import division, print_function, unicode_literals, absolute_import

import json
import os
import shutil
import unittest
import filecmp
from pymongo import MongoClient

from pymatgen.io.lammps.data import LammpsForceFieldData

from fireworks import LaunchPad, FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.lammps.workflows.core import get_wf_basic
from atomate.utils.testing import AtomateTest

__author__ = 'Kiran Mathew, Brandon Wood'
__email__ = 'kmathew@lbl.gov, b.wood@berkeley.edu'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "common", "test_files")

DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test
LAMMPS_CMD = None #"mpirun -n 4 lmp_mpi"


class TestLammpsWorkflows(AtomateTest):

    def setUp(self):
        super(TestLammpsWorkflows, self).setUp()
        self.data_file = os.path.join(module_dir, "test_files/peo.data")
        self.input_file = os.path.join(module_dir, "test_files/peo.in")
        self.db_file = os.path.join(db_dir, "db.json")
        self.reference_files_path = os.path.abspath(os.path.join(module_dir, "reference_files"))

    def test_lammps_wflow(self):

        if not LAMMPS_CMD:
            # fake run just copy files from reference path
            lammps_cmd = "cp -r {}/ .".format(self.reference_files_path)
        else:
            lammps_cmd = LAMMPS_CMD + " -in " + "lammps.in"
        print(lammps_cmd)
        wf = get_wf_basic("peo_test", self.input_file, self.data_file, "lammps.data", lammps_cmd,
                          is_forcefield=True, db_file=self.db_file)

        self.lp.add_wf(wf)
        # run
        rapidfire(self.lp, fworker=FWorker(env={"db_file": self.db_file}))
        d = self.get_task_collection().find_one()
        self._check_run(d)

    def _check_run(self, d):
        if not LAMMPS_CMD:
            path = d["dir_name"].split(":")[-1]
            self.assertTrue(filecmp.cmp(os.path.join(path, "lammps.log"),
                                        "{}/lammps.log".format(self.reference_files_path)))
            self.assertTrue(filecmp.cmp(os.path.join(path, "peo.dump"),
                                        "{}/peo.dump".format(self.reference_files_path)))
            self.assertTrue(filecmp.cmp(os.path.join(path, "peo.dcd"),
                                        "{}/peo.dcd".format(self.reference_files_path)))

if __name__ == "__main__":
    unittest.main()
