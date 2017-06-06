# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import shutil
import unittest

from pymatgen import Structure

from atomate.feff.fireworks.core import EXAFSPathsFW
from atomate.utils.testing import AtomateTest

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")


class TestFireworks(AtomateTest):

    def setUp(self):
        super(TestFireworks, self).setUp()
        self.struct = Structure.from_file(os.path.join(module_dir, "..", "..", "test_files",
                                                       "feo_781777.json"))

    def test_exafs_paths_fw(self):
        fw = EXAFSPathsFW(0, self.struct, [[249 , 0], [85, 0]])
        fw_dict = fw.as_dict()
        self.assertEqual(len(fw_dict["spec"]["_tasks"]), 5)
        ans = ['{{atomate.feff.firetasks.glue_tasks.CopyFeffOutputs}}',
               '{{atomate.feff.firetasks.write_inputs.WriteFeffFromIOSet}}',
               '{{atomate.feff.firetasks.write_inputs.WriteEXAFSPaths}}',
               '{{atomate.feff.firetasks.run_calc.RunFeffDirect}}',
               '{{atomate.feff.firetasks.parse_outputs.AddPathsToFilepadTask}}']
        self.assertEqual(ans, [ft["_fw_name"] for ft in fw_dict["spec"]["_tasks"]])


if __name__ == "__main__":
    unittest.main()
