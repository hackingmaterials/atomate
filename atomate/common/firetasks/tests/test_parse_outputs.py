
from __future__ import absolute_import, division, print_function, unicode_literals

import os
import unittest
import shutil

from fireworks import LaunchPad, Firework, Workflow
from fireworks.core.rocket_launcher import rapidfire
from pymatgen.apps.borg.hive import AbstractDrone
from atomate.common.firetasks.parse_outputs import ToDbTask

__author__ = 'Shyam Dwaraknath <shyamd@lbl.gov>'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "test_files")


class TestDrone(AbstractDrone):

    def assimilate(self, path):
        return {"Drone": "Test Drone"}

    def get_valid_paths(self, path):
        return path


class TestToDbTask(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.scratch_dir = os.path.join(module_dir, "scratch")

    def setUp(self):
        if os.path.exists(self.scratch_dir):
            shutil.rmtree(self.scratch_dir)
        os.makedirs(self.scratch_dir)
        os.chdir(self.scratch_dir)
        try:
            self.lp = LaunchPad.from_file(
                os.path.join(db_dir, "my_launchpad.yaml"))
            self.lp.reset("", require_password=False)

        except:
            raise unittest.SkipTest(
                'Cannot connect to MongoDB! Is the database server running? '
                'Are the credentials correct?')

    def tearDown(self):
        shutil.rmtree(self.scratch_dir)
        os.chdir(module_dir)
        self.lp.reset("", require_password=False)
        #self.lp.db.tasks.drop()

    def test_ToDbTask(self):
        # TODO: @shyamd Not sure how this test passed before? -computron
        raise unittest.SkipTest(
            'Shyam: please revise this test. The ToDbTask requires an actual drone=Drone not drone=str?')
        fw1 = Firework([ToDbTask(drone="atomate.common.firetasks.tests.test_parse_outputs.testDrone",
                                 mmdb="",
                                 db_file="db.json",
                                 calc_dir=db_dir)], name="fw1")

        wf = Workflow([fw1])
        self.lp.add_wf(wf)
        rapidfire(self.lp)

        task1 = self.lp.db.tasks.find_one({"task_id":1})
        self.assertEqual(task1['task_id'],1)
        self.assertEqual(task1['Drone'],'Test Drone')

if __name__ == "__main__":
    unittest.main()
