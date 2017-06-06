# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

from atomate.vasp.firetasks.glue_tasks import CopyVaspOutputs
from atomate.utils.testing import AtomateTest

__author__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

DEBUG_MODE = False


class TestCopyVaspOutputs(AtomateTest):

    @classmethod
    def setUpClass(cls):
        cls.plain_outdir = os.path.join(module_dir, "..", "..", "test_files",
                                        "Si_structure_optimization_plain", "outputs")
        cls.gzip_outdir = os.path.join(module_dir, "..", "..", "test_files",
                                       "Si_structure_optimization", "outputs")
        cls.relax2_outdir = os.path.join(module_dir, "..", "..", "test_files",
                                         "Si_structure_optimization_relax2", "outputs")

    def test_unittestsetup(self):
        files = ["INCAR", "KPOINTS", "POTCAR", "POSCAR", "CONTCAR", "OUTCAR"]
        for f in files:
            self.assertTrue(os.path.exists(os.path.join(self.plain_outdir, f)))
            self.assertTrue(os.path.exists(os.path.join(self.gzip_outdir, f + ".gz")))
            if f == "POTCAR":
                self.assertTrue(os.path.exists(os.path.join(self.relax2_outdir, f + ".gz")))
                self.assertTrue(os.path.exists(os.path.join(self.relax2_outdir, f + ".orig.gz")))
            else:
                self.assertTrue(os.path.exists(os.path.join(self.relax2_outdir, f + ".relax1.gz")))
                self.assertTrue(os.path.exists(os.path.join(self.relax2_outdir, f + ".relax2.gz")))

    def test_plain_copy(self):
        ct = CopyVaspOutputs(calc_dir=self.plain_outdir)
        ct.run_task({})
        files = ["INCAR", "KPOINTS", "POSCAR", "POTCAR", "OUTCAR"]
        for f in files:
            self.assertTrue(os.path.exists(os.path.join(self.scratch_dir, f)))

        no_files = ["CONTCAR", "EIGENVAL", "IBZKPT"]
        for f in no_files:
            self.assertFalse(os.path.exists(os.path.join(self.scratch_dir, f)))

        # make sure CONTCAR was copied properly
        with open(os.path.join(self.plain_outdir, "CONTCAR")) as f1:
            with open(os.path.join(self.scratch_dir, "POSCAR")) as f2:
                self.assertEqual(f1.read(), f2.read())

    def test_plain_copy_more(self):
        ct = CopyVaspOutputs(calc_dir=self.plain_outdir,
                             contcar_to_poscar=False,
                             additional_files=["IBZKPT"])
        ct.run_task({})
        files = ["INCAR", "KPOINTS", "POSCAR", "POTCAR", "IBZKPT"]
        for f in files:
            self.assertTrue(os.path.exists(os.path.join(self.scratch_dir, f)))

        no_files = ["CONTCAR", "EIGENVAL"]
        for f in no_files:
            self.assertFalse(os.path.exists(os.path.join(self.scratch_dir, f)))

        # make sure CONTCAR was NOT copied and POSCAR was instead copied
        with open(os.path.join(self.plain_outdir, "POSCAR")) as f1:
            with open(os.path.join(self.scratch_dir, "POSCAR")) as f2:
                self.assertEqual(f1.read(), f2.read())

    def test_gzip_copy(self):
        ct = CopyVaspOutputs(calc_dir=self.gzip_outdir)
        ct.run_task({})
        files = ["INCAR", "KPOINTS", "POTCAR", "POSCAR"]
        for f in files:
            self.assertTrue(os.path.exists(os.path.join(self.scratch_dir, f)))

        no_files = ["CONTCAR", "EIGENVAL"]
        for f in no_files:
            self.assertFalse(os.path.exists(os.path.join(self.scratch_dir, f)))

    def test_relax2_copy(self):
        ct = CopyVaspOutputs(calc_dir=self.relax2_outdir, additional_files=["IBZKPT"])
        ct.run_task({})
        files = ["INCAR", "KPOINTS", "POTCAR", "POSCAR", "IBZKPT"]
        for f in files:
            self.assertTrue(os.path.exists(os.path.join(self.scratch_dir, f)))

        no_files = ["CONTCAR", "EIGENVAL"]
        for f in no_files:
            self.assertFalse(os.path.exists(os.path.join(self.scratch_dir, f)))


if __name__ == "__main__":
    unittest.main()
