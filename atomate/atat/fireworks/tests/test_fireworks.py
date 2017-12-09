# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import shutil
import unittest

from pymatgen import Structure

from atomate.atat.fireworks.core import McsqsFW
from atomate.utils.testing import AtomateTest

__author__ = 'Matthew Horton'
__email__ = 'mkhorton@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")


class TestFireworks(AtomateTest):

    def setUp(self):
        super(TestFireworks, self).setUp()
        self.struct = Structure.from_file(os.path.join(module_dir, "..", "..", "test_files",
                                                       "MnSrCaO3.json"))

    def test_mcsqs_fw(self):
        fw = McsqsFW(self.struct)
        fw_dict = fw.as_dict()
        self.assertEqual(fw_dict['spec']['_tasks'][0]['_fw_name'], 'ScriptTask')


if __name__ == "__main__":
    unittest.main()
