# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import os
import shutil
import unittest

from matmethods.vasp.firetasks.glue_tasks import CopyVaspOutputs

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
scratch_dir = os.path.join(module_dir, "scratch")

DEBUG_MODE = False


class TestCopyVaspOutputs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.plain_outdir = os.path.join(module_dir, "reference_files",
                                        "Si_structure_optimization", "outputs")
        cls.gzip_outdir = os.path.join(module_dir, "reference_files",
                                       "Si_nscf_uniform", "outputs")

    def setUp(self):
        if os.path.exists(scratch_dir):
            shutil.rmtree(scratch_dir)
        os.makedirs(scratch_dir)
        os.chdir(scratch_dir)

    def tearDown(self):
        if not DEBUG_MODE:
            shutil.rmtree(scratch_dir)

    def test_unittestsetup(self):
        files = ["INCAR", "KPOINTS", "POTCAR", "POSCAR", "CONTCAR", "OUTCAR"]
        for f in files:
            self.assertTrue(os.path.exists(os.path.join(self.plain_outdir, f)))
            self.assertTrue(
                os.path.exists(os.path.join(self.gzip_outdir, f + ".gz")))

    def test_plain_copy(self):
        ct = CopyVaspOutputs(vasp_dir=self.plain_outdir)
        ct.run_task({})
        files = ["INCAR", "KPOINTS", "POSCAR", "POTCAR", "OUTCAR"]
        for f in files:
            self.assertTrue(os.path.exists(os.path.join(scratch_dir, f)))

        no_files = ["CONTCAR", "EIGENVAL", "IBZKPT"]
        for f in no_files:
            self.assertFalse(os.path.exists(os.path.join(scratch_dir, f)))

        # make sure CONTCAR was copied properly
        with open(os.path.join(self.plain_outdir, "CONTCAR")) as f1:
            with open(os.path.join(scratch_dir, "POSCAR")) as f2:
                self.assertEqual(f1.read(), f2.read())

    def test_plain_copy_more(self):
        ct = CopyVaspOutputs(vasp_dir=self.plain_outdir,
                             contcar_to_poscar=False,
                             additional_files=["IBZKPT"])
        ct.run_task({})
        files = ["INCAR", "KPOINTS", "POSCAR", "POTCAR", "IBZKPT"]
        for f in files:
            self.assertTrue(os.path.exists(os.path.join(scratch_dir, f)))

        no_files = ["CONTCAR", "EIGENVAL"]
        for f in no_files:
            self.assertFalse(os.path.exists(os.path.join(scratch_dir, f)))

        # make sure CONTCAR was NOT copied and POSCAR was instead copied
        with open(os.path.join(self.plain_outdir, "POSCAR")) as f1:
            with open(os.path.join(scratch_dir, "POSCAR")) as f2:
                self.assertEqual(f1.read(), f2.read())

    def test_gzip_copy(self):
        ct = CopyVaspOutputs(vasp_dir=self.gzip_outdir)
        ct.run_task({})
        files = ["INCAR", "KPOINTS", "POTCAR", "POSCAR"]
        for f in files:
            self.assertTrue(os.path.exists(os.path.join(scratch_dir, f)))

        no_files = ["CONTCAR", "EIGENVAL"]
        for f in no_files:
            self.assertFalse(os.path.exists(os.path.join(scratch_dir, f)))


if __name__ == "__main__":
    unittest.main()
