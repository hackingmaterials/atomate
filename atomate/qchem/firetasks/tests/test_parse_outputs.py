# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
import shutil

from atomate.qchem.firetasks.parse_outputs import QChemToDb
from atomate.utils.testing import AtomateTest

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "11/14/18"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")


class TestParseOutputQChem(AtomateTest):

    def setUp(self, lpad=False):
        super(TestParseOutputQChem, self).setUp(lpad=False)

    def tearDown(self):
        pass

    def test_parse_ion_placing_FF(self):
        with patch("atomate.qchem.firetasks.parse_outputs.FWAction"
                   ) as FWAction_patch:
            ft = QChemToDb(
                    calc_dir=os.path.join(module_dir, "..", "..", "test_files", "ion_placer_files", "launcher0"),
                    db_file=os.path.join(db_dir, "db.json"),
                    input_file="mol.qin",
                    output_file="mol.qout",
                    additional_fields={
                        "task_label": "ion_pos_0",
                        "special_run_type": "frequency_flattener",
                        "linked": False
                    })
            ft.run_task({})
            self.assertEqual(FWAction_patch.call_args[1]["update_spec"]["ion_pos_0"]["energy"],-349.91898078154)


if __name__ == "__main__":
    unittest.main()
