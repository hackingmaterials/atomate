# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure

from atomate.vasp.fireworks.nmr import *

from atomate.vasp.firetasks.write_inputs import WriteVaspNMRFromPrev, WriteVaspFromIOSet

__author__ = 'Shyam Dwaraknath'
__email__ = 'shyamd@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
reference_dir = os.path.join(module_dir, "..", "..", "test_files")

DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...


class TestNMRFireworks(unittest.TestCase):
    def setUp(self):
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = Lattice([[3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        self.structure = Structure(lattice, ["Si", "Si"], coords)

    def test_NMRFW(self):
        self.assertEqual(NMRFW(structure=self.structure).name, "Si-nmr tensor")

        # Ensure correct Fireworks structure based on input type
        # Simple structure input
        self.assertEqual(type(NMRFW(structure=self.structure).tasks[0]), WriteVaspFromIOSet)

        # Previous calc dir input
        self.assertEqual(type(NMRFW(prev_calc_dir="/some_random_dir").tasks[1]), WriteVaspNMRFromPrev)
        self.assertEqual(NMRFW(prev_calc_dir="/some_random_dir").tasks[0]["calc_dir"], "/some_random_dir")

        # Direct to parent calculations
        self.assertEqual(type(NMRFW(copy_vasp_outputs=True, parents=[None]).tasks[1]), WriteVaspNMRFromPrev)
        self.assertEqual(NMRFW(copy_vasp_outputs=True, parents=[None]).tasks[0]["calc_loc"], True)

        self.assertRaises(ValueError, NMRFW)


if __name__ == "__main__":
    unittest.main()
