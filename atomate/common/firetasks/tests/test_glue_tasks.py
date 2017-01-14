from __future__ import absolute_import

import os
import shutil
import unittest

from fireworks import LaunchPad
from fireworks.core.firework import Firework, Workflow
from fireworks.core.rocket_launcher import rapidfire
from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.utils.utils import get_calc_loc

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "reference_files",
                      "db_connections")


class TestPassCalcLocs(unittest.TestCase):

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
        self.lp.reset("", require_password=False)

    def test_passcalclocs(self):
        fw1 = Firework([PassCalcLocs(name="fw1")], name="fw1")
        fw2 = Firework([PassCalcLocs(name="fw2")], name="fw2", parents=fw1)
        fw3 = Firework([PassCalcLocs(name="fw3")], name="fw3", parents=fw2)

        wf = Workflow([fw1, fw2, fw3])
        self.lp.add_wf(wf)
        rapidfire(self.lp)

        fw2 = self.lp.get_fw_by_id(self.lp.get_fw_ids({"name": "fw2"})[0])
        fw3 = self.lp.get_fw_by_id(self.lp.get_fw_ids({"name": "fw3"})[0])

        self.assertEqual(len(fw2.spec["calc_locs"]), 1)
        self.assertEqual(len(fw3.spec["calc_locs"]), 2)
        self.assertEqual(fw3.spec["calc_locs"][0]["name"], "fw1")
        self.assertEqual(fw3.spec["calc_locs"][1]["name"], "fw2")
        self.assertNotEqual(fw3.spec["calc_locs"][0]["path"],
                            fw3.spec["calc_locs"][1]["path"])

        calc_locs = fw3.spec["calc_locs"]
        self.assertEqual(get_calc_loc("fw1", calc_locs), calc_locs[0])
        self.assertEqual(get_calc_loc("fw2", calc_locs), calc_locs[1])
        self.assertEqual(get_calc_loc(True, calc_locs), calc_locs[1])



