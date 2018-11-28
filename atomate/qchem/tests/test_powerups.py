# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
import shutil

from atomate.qchem.firetasks.write_inputs import WriteInputFromIOSet
from atomate.qchem.firetasks.run_calc import RunQChemDirect
from atomate.qchem.firetasks.parse_outputs import QChemToDb
from fireworks import Firework, Workflow
from atomate.utils.testing import AtomateTest
from pymatgen.io.qchem.outputs import QCOutput
from atomate.qchem.powerups import use_fake_qchem

__author__ = "Brandon Wood"
__email__ = "b.wood@berkeley.edu"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")


class TestPowerups(AtomateTest):
    @classmethod
    def setUpClass(cls):
        out_file = os.path.join(module_dir, "..", "test_files", "FF_working",
                                "test.qout.opt_1")
        qc_out = QCOutput(filename=out_file)
        cls.act_mol = qc_out.data["molecule_from_optimized_geometry"]

    def setUp(self, lpad=False):
        super(TestPowerups, self).setUp(lpad=False)

    def tearDown(self):
        # this removes the scratch dir made by AtomateTest
        shutil.rmtree(self.scratch_dir)
        # this removes the file that gets written
        for x in ["mol.qin", "task.json"]:
            if os.path.exists(os.path.join(module_dir, x)):
                os.remove(os.path.join(module_dir, x))

    def test_use_fake_qchem(self):

        input_file = "test.qin.opt_1"
        output_file = "test.qout.opt_1"
        calc_dir = os.path.join(module_dir, "..", "test_files", "FF_working")

        run_task = RunQChemDirect(
            qchem_cmd="should not need this going to be replaced with fake run"
        )
        p_task = QChemToDb(
            calc_dir=calc_dir, input_file=input_file, output_file=output_file)
        fw1 = Firework([run_task, p_task], name="test_fake_run")
        w_task = WriteInputFromIOSet(
            qchem_input_set="OptSet", write_to_dir=module_dir)
        fw2 = Firework([w_task], parents=fw1, name="test_write")
        wf = Workflow([fw1, fw2])
        ref_dirs = {"test_fake_run": os.path.join(calc_dir, output_file)}

        fake_run_wf = use_fake_qchem(wf, ref_dirs)
        test_fake_run = False
        for fw in fake_run_wf.fws:
            if fw.name == "test_fake_run":
                for t in fw.tasks:
                    if "RunQChemFake" in str(t):
                        test_fake_run = True
        self.assertTrue(test_fake_run)


if __name__ == "__main__":
    unittest.main()
