# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

from atomate.qchem.firetasks.write_inputs import WriteInputFromIOSet
from atomate.utils.testing import AtomateTest
from pymatgen.core import Molecule
from pymatgen.io.qchem_io.inputs import QCInput

__author__ = 'Brandon Wood'
__email__ = 'b.wood@berkeley.edu'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class TestWriteInputQChem(AtomateTest):
    @classmethod
    def setUpClass(cls):

        co_species = ['C', 'O']
        co_coords = [[0.0, 0.0, 0.0], [1.3, 0.0, 0.0]]
        cls.co_mol = Molecule(co_species, co_coords)
        cls.co_opt_ref_in = QCInput.from_file(
            os.path.join(module_dir, "..", "..", "test_files", "co_qc.in"))

    def setUp(self, lpad=False):
        super(TestWriteInputQChem, self).setUp(lpad=False)

    def tearDown(self):
        for x in ["qc.in"]:
            if os.path.exists(os.path.join(module_dir, x)):
                os.remove(os.path.join(module_dir, x))

    def test_write_input_from_io_set(self):
        ft = WriteInputFromIOSet({"molecule": self.co_mol, "qchem_input_set": "OptSet"})
        ft.run_task({})
        self.assertEqual(str(QCInput.from_file("qc.in")), str(self.co_opt_ref_in))

if __name__ == '__main__':
    unittest.main()
