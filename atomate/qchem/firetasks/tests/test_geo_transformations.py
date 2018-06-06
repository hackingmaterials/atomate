# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

from atomate.qchem.firetasks.geo_transformations import RotateTorsion
from atomate.utils.testing import AtomateTest
from pymatgen.core import Molecule
import numpy as np

__author__ = "Brandon Wood"
__email__ = "b.wood@berkeley.edu"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class TestGeoTransformations(AtomateTest):
    @classmethod
    def setUpClass(cls):

        cls.pt_mol = Molecule.from_file(
            os.path.join(module_dir, "..", "..", "test_files",
                         "pt_gs_wb97mv_tz_initial.xyz"))
        cls.pt_rot_90_mol = Molecule.from_file(
            os.path.join(module_dir, "..", "..", "test_files",
                         "pt_rotated_90.0.xyz"))

    def setUp(self, lpad=False):
        super(TestGeoTransformations, self).setUp(lpad=False)

    def tearDown(self):
        pass

    def test_rotate_torsion(self):
        atom_indexes = [6, 8, 9, 10]
        angle = 90.0
        ft = RotateTorsion({
            "molecule": self.pt_mol,
            "atom_indexes": atom_indexes,
            "angle": angle
        })
        rot_mol = ft.run_task({})
        test_mol = Molecule.from_dict(
            rot_mol.as_dict()["update_spec"]["prev_calc_molecule"])
        np.testing.assert_equal(self.pt_rot_90_mol.species, test_mol.species)
        np.testing.assert_allclose(
            self.pt_rot_90_mol.cart_coords, test_mol.cart_coords, atol=0.0001)


if __name__ == "__main__":
    unittest.main()
