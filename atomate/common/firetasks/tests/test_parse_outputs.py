
from __future__ import absolute_import, division, print_function, unicode_literals

import os
import unittest
import shutil

from fireworks import LaunchPad, Firework, Workflow
from fireworks.core.rocket_launcher import rapidfire

from pymatgen.apps.borg.hive import AbstractDrone

from atomate.common.firetasks.parse_outputs import ToDbTask
from atomate.utils.testing import AtomateTest

__author__ = 'Shyam Dwaraknath <shyamd@lbl.gov>'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "test_files")


class TestDrone(AbstractDrone):

    def __init__(self):
        pass

    def assimilate(self, path):
        return {"drone": "Test Drone",
                "dir_name": "/test",
                'state': "successful"}

    def get_valid_paths(self, path):
        return path


class TestToDbTask(AtomateTest):

    def test_ToDbTask(self):
        raise unittest.SkipTest("Shyam - please fix this test.")
        d = TestDrone()

        fw1 = Firework([ToDbTask(drone=d,
                                 mmdb="atomate.vasp.database.VaspCalcDb",
                                 db_file=os.path.join(db_dir, "db.json"),
                                 calc_dir=db_dir)], name="fw1")

        wf = Workflow([fw1])
        self.lp.add_wf(wf)
        rapidfire(self.lp)

        task1 = self.lp.db.tasks.find_one({"task_id":1})
        self.assertEqual(task1['task_id'],1)
        self.assertEqual(task1['dir_name'],"/test" )
        self.assertEqual(task1['drone'],'Test Drone')

if __name__ == "__main__":
    unittest.main()
