# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import shutil
import unittest
from glob import glob

from pymatgen import Structure
from pymatgen.io.feff.sets import MPEXAFSSet
from pymatgen.io.feff.inputs import Paths

from atomate.feff.firetasks.glue_tasks import CopyFeffOutputs
from atomate.feff.firetasks.write_inputs import WriteEXAFSPaths

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")


class TestTasks(unittest.TestCase):

    def setUp(self):
        self.struct = Structure.from_file(os.path.join(module_dir, "..", "..", "test_files",
                                                       "feo_781777.json"))
        self.scratch_dir = os.path.join(module_dir, "scratch")
        if os.path.exists(self.scratch_dir):
            shutil.rmtree(self.scratch_dir)
        os.makedirs(self.scratch_dir)
        os.chdir(self.scratch_dir)

    def test_copy_feff_outputs_task(self):
        t = CopyFeffOutputs(calc_dir=os.path.join(module_dir, "..", "..", "test_files"))
        t.run_task({})
        ans = ["Co2O2.cif", "feff_eels.inp", "feo_781777.json"]
        files = glob("*")
        self.assertEqual(sorted(ans), sorted(files))

    def test_write_paths_task(self):
        exafs = MPEXAFSSet(0, self.struct, edge='K', radius=10)
        t = WriteEXAFSPaths(feff_input_set=exafs, paths=[[249 , 0], [85, 0]])
        paths = Paths(exafs.atoms, [[249, 0], [85, 0]])
        paths.write_file("paths_ans.dat")
        t.run_task({})
        with open("paths_ans.dat", "r") as ans, open("paths.dat", "r") as tmp:
            self.assertEqual(ans.readlines(), tmp.readlines())

    def tearDown(self):
        shutil.rmtree(self.scratch_dir)


if __name__ == "__main__":
    unittest.main()
