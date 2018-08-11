# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

from fireworks import FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.vasp.powerups import use_fake_vasp
from atomate.vasp.workflows.base.grain_boundary import get_wf_gb
from atomate.utils.testing import AtomateTest
from pymatgen.core import Structure
from pymatgen.analysis.gb.gb import GBGenerator

__author__ = 'Hui Zheng'
__email__ = 'huz071@eng.ucsd.edu'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
ref_dir = os.path.join(module_dir, "..", "..", "test_files")

DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...

gb_gen_params_s5 = {"rotation_axis": [1, 0, 0], "rotation_angle": 36.86989764584402,
                    "expand_times": 1, "vacuum_thickness": 0.0, "normal": True,
                    "ratio": None, "plane": None}

gb_gen_params_s3 = {"rotation_axis": [1, 1, 1], "rotation_angle": 60.0,
                    "expand_times": 2, "vacuum_thickness": 0.0, "normal": True,
                    "ratio": None, "plane": None}


class TestGrainboundaryWorkflow(AtomateTest):
    def setUp(self, lpad=True):
        super(TestGrainboundaryWorkflow, self).setUp()
        self.bulk_structure = Structure.from_spacegroup("Im-3m", Lattice.cubic(3.42682), ["Li", "Li"],
                                                        [[0, 0, 0], [0.5, 0.5, 0.5]])
        gbg = GBGenerator(self.bulk_structure)
        gb_100_s5 = gbg.gb_from_parameters(**gb_gen_params_s5)
        gb_info_s5 = gb_100_s5.as_dict()
        self.wf_1 = get_wf_gb(gb=None,
                              bulk_structure=self.bulk_structure, gb_gen_params=gb_gen_params_s3,
                              db_file=os.path.join(db_dir, "db.json"))
        self.wf_2 = get_wf_gb(gb=gb_100_s5, bulk_structure=None, gb_gen_params=None,
                              additional_info=gb_info_s5, tag=["mp-135"],
                              db_file=os.path.join(db_dir, "db.json"))

    @classmethod
    def _simulate_vasprun_start_from_bulk(cls, wf):
        reference_dir = os.path.abspath(os.path.join(ref_dir, "grain_boundary_wf"))
        li_ref_dir = {"Li-bulk-structure optimization": os.path.join(reference_dir, "1"),
                      "Li_gb_s3_(1, 1, 1) optimization": os.path.join(reference_dir, "2")}
        return use_fake_vasp(wf, li_ref_dir, params_to_check=["ENCUT", "ISIF", "IBRION"])

    def _check_run(self, d, mode):
        if mode not in ["Li_gb_s3_(1, 1, 1) optimization", "Li_gb_s5_(1, 0, 0) optimization",
                        "Li-bulk-structure optimization"]:
            raise ValueError("Invalid mode!")
        if "optimization" in mode:
            self.assertEqual(d["formula_reduced_abc"], "Li1")

    def test_start_from_bulk_wf(self):
        wf = self._simulate_vasprun_start_from_bulk(self.wf_1)
        self.assertEqual(len(self.wf_1.fws), 2)
        self.lp.add_wf(wf)
        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))

        # check relaxation
        d = self.get_task_collection().find_one({"task_label": "Li_gb_s3_(1, 1, 1) optimization"})
        self._check_run(d, mode="Li_gb_s3_(1, 1, 1) optimization")

        wf = self.lp.get_wf_by_fw_id(1)
        self.assertTrue(all([s == 'COMPLETED' for s in wf.fw_states.values()]))

    @classmethod
    def _simulate_vasprun_start_from_gb(cls, wf):
        reference_dir = os.path.abspath(os.path.join(ref_dir, "grain_boundary_wf"))
        li_ref_dir = {"Li_gb_s5_(1, 0, 0) optimization": os.path.join(reference_dir, "3")}
        return use_fake_vasp(wf, li_ref_dir, params_to_check=["ENCUT", "ISIF", "IBRION"])

    def test_start_from_gb_wf(self):
        wf = self._simulate_vasprun_start_from_gb(self.wf_2)
        self.assertEqual(len(self.wf_1.fws), 1)
        self.lp.add_wf(wf)
        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))

        # check relaxation
        d = self.get_task_collection().find_one({"task_label": "Li_gb_s5_(1, 0, 0) optimization"})
        print(d.keys(), "********test what is the keys of d *****")
        self._check_run(d, mode="Li_gb_s5_(1, 0, 0) optimization")

        wf = self.lp.get_wf_by_fw_id(0)
        self.assertTrue(all([s == 'COMPLETED' for s in wf.fw_states.values()]))


if __name__ == "__main__":
    unittest.main()
