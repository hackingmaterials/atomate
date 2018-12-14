# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

import numpy as np

from fireworks import FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.vasp.powerups import use_fake_vasp
from atomate.vasp.workflows.presets.core import wf_nmr
from atomate.utils.testing import AtomateTest

from pymatgen import Structure

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
reference_dir = os.path.join(module_dir, "..", "..", "test_files", "nmr_wf")

DEBUG_MODE = False  # If True, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...
_write_task_docs = False  # Test developer option: defaults to False, need to be True only once


class TestNMRWorkflow(AtomateTest):
    """
    This test will either actually run VASP (if VASP_CMD is set) or artificially pass on outputs
    (if not VASP_CMD) and test the NMR workflow and its implementation and outputs
    for an example calculation for LiAlSiO4.
    """

    def setUp(self):
        super(TestNMRWorkflow, self).setUp()
        self.struct = Structure.from_file(os.path.join(reference_dir, "LiAlSiO4.json"))
        self.wf = wf_nmr(self.struct)

    def test_wf(self):
        self.wf = self._simulate_vasprun(self.wf)

        self.assertEqual(len(self.wf.fws), 3)

        self.lp.add_wf(self.wf)

        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))

        # check relaxation
        d = self.get_task_collection().find_one({"task_label": "structure optimization"})
        self._check_run(d, mode="structure optimization")
        # check two of the deformation calculations
        d = self.get_task_collection().find_one({"task_label": "cs tensor"})
        self._check_run(d, mode="cs tensor")

        d = self.get_task_collection().find_one({"task_label": "efg tensor"})
        self._check_run(d, mode="efg tensor")

        wf = self.lp.get_wf_by_fw_id(1)
        self.assertTrue(all([s == 'COMPLETED' for s in wf.fw_states.values()]))

    def _simulate_vasprun(self, wf):
        ref_dirs = {
            "structure optimization": os.path.join(reference_dir, "structure_optimization"),
            "cs tensor": os.path.join(reference_dir, "cs"),
            "efg tensor": os.path.join(reference_dir, "efg")
        }

        if not VASP_CMD:
            wf = use_fake_vasp(wf, ref_dirs)
        else:
            wf = use_custodian(wf)

        return wf

    def _check_run(self, d, mode):
        if mode not in ["structure optimization", "cs tensor", "efg tensor"]:
            raise ValueError("Invalid mode!")

        self.assertEqual(d["formula_pretty"], "LiAlSiO4")
        self.assertEqual(d["formula_anonymous"], "ABCD4")
        self.assertEqual(d["nelements"], 4)
        self.assertEqual(d["state"], "successful")

        if mode == "structure optimization":
            self.assertAlmostEqual(d["calcs_reversed"][0]["output"]["structure"]["lattice"]["a"], 5.297, 2)
            self.assertAlmostEqual(d["output"]["energy_per_atom"], -7.226, 2)

        if mode == "cs tensor":
            cs = d["calcs_reversed"][0]["output"]["outcar"]["nmr_cs"]
            np.testing.assert_allclose(
                cs["g0"],
                [[-3.896026, 0.008131, 0.004022], [0.010171, -4.303434, 0.007171], [0.013091, 0.007672, -4.478203]])
            np.testing.assert_allclose(cs["raw"][0],
                                       [[-56.972805, -0.000285, -0.000021], [-0.000354, -47.358647, 0.000055],
                                        [-0.000065, 0.000144, -47.036531]])

            np.testing.assert_allclose(cs["raw"][-1],
                                       [[-13.001104, -20.568634, 19.91829], [-25.384874, -37.953075, 28.781717],
                                        [24.920975, 29.464446, -29.6856]])
        elif mode == "efg tensor":
            efg = d["calcs_reversed"][0]["output"]["outcar"]["nmr_efg"]
            np.testing.assert_allclose(efg["raw"][0],
                                       [[-0.428, -0.105, -0.0], [-0.105, 0.388, -0.0], [-0.0, -0.0, 0.041]])
            np.testing.assert_allclose(
                efg["raw"][-1], [[3.145, 29.549, -26.609], [29.549, 4.044, -38.683], [-26.609, -38.683, -7.188]])


if __name__ == "__main__":
    unittest.main()
