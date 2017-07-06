# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import json
import os
import unittest

import numpy as np
from monty.json import MontyEncoder

from fireworks import FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.vasp.powerups import use_fake_vasp, use_no_vasp
from atomate.vasp.workflows.presets.core import wf_bulk_modulus
from atomate.utils.testing import AtomateTest

from pymatgen import Structure
from pymatgen.util.testing import PymatgenTest
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
reference_dir = os.path.join(module_dir, "..", "..", "test_files", "bulk_modulus_wf")

DEBUG_MODE = False  # If True, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...
_write_task_docs = False # Test developer option: defaults to False, need to be True only once


class TestBulkModulusWorkflow(AtomateTest):
    """
    This test will either actually run VASP (if VASP_CMD is set) or artificially pass on outputs
    (if not VASP_CMD) and test the whole bulk modulus workflow and its implementation and outputs
    for an example calculation for silicon.

    note for the developer of the test:
    This tests can be run in two modes if not VASP_CMD:
    1. first all inputs and outputs of all deformations are present in which case
        in each folder a trimmed version of task.json will be generated so that
    2. once task.json is present, VaspRun can be skipped where task.json is
        available and their "inputs" and "outputs" folders can be removed
    """

    def setUp(self):
        super(TestBulkModulusWorkflow, self).setUp()
        self.struct_si = PymatgenTest.get_structure("Si")
        self.ndeformations = 6
        self.deformations = [(np.identity(3) * (1 + x)).tolist() for x in
                             np.linspace(-0.05, 0.05, self.ndeformations)]
        self.wf_config = {"VASP_CMD": ">>vasp_cmd<<", "DB_FILE": ">>db_file<<"}
        self.wf = wf_bulk_modulus(self.struct_si, self.wf_config)

    def _simulate_vasprun(self, wf):
        no_vasp_ref_dirs = {}
        fake_vasp_ref_dirs = {}
        for i in range(2, self.ndeformations+2):
            if os.path.exists(os.path.join(reference_dir, str(i), "inputs")):
                if not VASP_CMD:
                    fake_vasp_ref_dirs["bulk_modulus deformation {}".format(i-2)] = os.path.join(reference_dir, str(i))
            else:
                no_vasp_ref_dirs["bulk_modulus deformation {}".format(i-2)] = os.path.join(self.scratch_dir, str(i))

        fake_vasp_ref_dirs["structure optimization"] = os.path.join(reference_dir, "1")
        new_wf = use_no_vasp(wf, no_vasp_ref_dirs)
        return use_fake_vasp(new_wf, fake_vasp_ref_dirs, params_to_check=["ENCUT"])

    def _check_run(self, d, mode):
        if mode not in ["structure optimization", "bulk_modulus deformation 0",
                        "bulk_modulus deformation 4", "fit equation of state"]:
            raise ValueError("Invalid mode!")

        if mode not in ["fit equation of state"]:
            self.assertEqual(d["formula_pretty"], "Si")
            self.assertEqual(d["formula_anonymous"], "A")
            self.assertEqual(d["nelements"], 1)
            self.assertEqual(d["state"], "successful")

        if mode in ["structure optimization"]:
            self.assertAlmostEqual(d["calcs_reversed"][0]["output"]["structure"]["lattice"]["a"], 3.866, 3)
            self.assertAlmostEqual(d["output"]["energy_per_atom"], -5.432, 3)
            self.relaxed_struct_si = d["calcs_reversed"][0]["output"]["structure"]

        elif mode in ["bulk_modulus deformation 0"]:
            for i, l in enumerate(["a", "b", "c"]):
                self.assertAlmostEqual(d["input"]["structure"]["lattice"][l],
                self.relaxed_struct_si["lattice"][l]* (self.deformations[0][i][i]), 2)
            stress = d["calcs_reversed"][0]["output"]["ionic_steps"][-1]["stress"]
            np.testing.assert_allclose(stress, np.diag([189.19, 189.19, 189.19]), atol=1e-2)

        elif mode in ["bulk_modulus deformation 4"]:
            for i, l in enumerate(["a", "b", "c"]):
                self.assertAlmostEqual(d["input"]["structure"]["lattice"][l],
                        self.relaxed_struct_si["lattice"][l] * (self.deformations[4][i][i]), 2)
            stress = d["calcs_reversed"][0]["output"]["ionic_steps"][-1]["stress"]
            np.testing.assert_allclose(stress, np.diag([-65.56, -65.56, -65.56]), atol=1e-2)

        elif mode in ["fit equation of state"]:
            self.assertAlmostEqual(d["bulk_modulus"], 88.90, places=2)
            self.assertEqual(len(d["all_task_ids"]), 7)
            self.assertEqual(len(d["energies"]), self.ndeformations)
            self.assertEqual(len(d["volumes"]), self.ndeformations)
            s = SpacegroupAnalyzer(Structure.from_dict(d["structure"])).get_conventional_standard_structure()
            self.assertAlmostEqual(s.lattice.c, 5.468, places=3)

    def setup_task_docs(self):
        self.task_file = "task.json"
        for i in range(2, self.ndeformations+2):
            if os.path.exists(os.path.join(reference_dir, str(i), self.task_file)):
                with open(os.path.join(reference_dir, str(i), self.task_file)) as fp:
                    d = json.load(fp)
                    new_fw = self.lp.fireworks.find_one(
                        {"name": {"$regex": "bulk_modulus deformation {}".format(i - 2)}})

                    # the fw tag (inluded in "name") is important in pulling tasks in the last FW
                    # in wf_bulk_modulus
                    d["task_label"] = new_fw["name"]
                    d["task_id"] += i + 1000000  # to avoid duplicate task_id

                os.makedirs(os.path.join(self.scratch_dir, str(i)))
                with open(os.path.join(self.scratch_dir, str(i), "task.json"), 'w') as fp:
                    json.dump(d, fp, sort_keys=True, indent=4, ensure_ascii=False, cls=MontyEncoder)

            elif not os.path.exists(os.path.join(reference_dir, str(i), "inputs")):
                raise IOError("neither {} nor {} are present in {}".format(
                    "inputs", self.task_file, os.path.join(reference_dir, str(i))))

    def write_task_docs(self):
        # this step needs to be run once: once task.json is present, remove the inputs/outputs folders
        for i in range(2, self.ndeformations + 2):
            # not to unnecessarily override available task.json
            if not os.path.exists(os.path.join(reference_dir, str(i), "task.json")):
                d = self.get_task_collection().find_one(
                    {"task_label": {"$regex": "bulk_modulus deformation {}".format(i-2)}})
                rm_props = ["bandstructure", "input"]
                for icalc in range(len(d["calcs_reversed"])):
                    for prop in rm_props:
                        try:
                            del (d["calcs_reversed"][icalc][prop])
                        except:
                            pass
                with open(os.path.join(reference_dir, str(i), "task.json"), 'w') as fp:
                    json.dump(d, fp, sort_keys=True, indent=4, ensure_ascii=False, cls=MontyEncoder)

    def test_wf(self):
        self.wf = self._simulate_vasprun(self.wf)
        self.assertEqual(len(self.wf.fws), self.ndeformations+2)

        defo_vis = [fw.tasks[2]['vasp_input_set'] for fw in self.wf.fws if "deform" in fw.name]
        assert all([vis.user_incar_settings['NSW'] == 99 for vis in defo_vis])
        assert all([vis.user_incar_settings['IBRION'] == 2 for vis in defo_vis])

        self.lp.add_wf(self.wf)

        # this is specific to bulk_modulus_wf "fit equation of state" that uses FW tag
        self.setup_task_docs()

        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))
        if _write_task_docs:
            self.write_task_docs()

        # check relaxation
        d = self.get_task_collection().find_one({"task_label": {"$regex": "structure optimization"}})
        self._check_run(d, mode="structure optimization")

        # check two of the deformation calculations
        d = self.get_task_collection().find_one({"task_label": {"$regex": "bulk_modulus deformation 0"}})
        self._check_run(d, mode="bulk_modulus deformation 0")

        d = self.get_task_collection().find_one({"task_label": {"$regex": "bulk_modulus deformation 4"}})
        self._check_run(d, mode="bulk_modulus deformation 4")

        # check the final results
        d = self.get_task_collection(coll_name="eos").find_one()
        self._check_run(d, mode="fit equation of state")


if __name__ == "__main__":
    unittest.main()
