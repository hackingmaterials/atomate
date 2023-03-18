import os
import unittest

import numpy as np
from pymatgen.core.structure import Molecule
from pymatgen.io.qchem.outputs import QCOutput

from atomate.qchem.firetasks.geo_transformations import PerturbGeometry, RotateTorsion
from atomate.utils.testing import AtomateTest

__author__ = "Brandon Wood, Evan Spotte-Smith"
__email__ = "b.wood@berkeley.edu"

module_dir = os.path.dirname(os.path.abspath(__file__))


class TestGeoTransformations(AtomateTest):
    @classmethod
    def setUpClass(cls):

        cls.pt_mol = Molecule.from_file(
            os.path.join(
                module_dir, "..", "..", "test_files", "pt_gs_wb97mv_tz_initial.xyz"
            )
        )
        cls.pt_rot_90_mol = Molecule.from_file(
            os.path.join(module_dir, "..", "..", "test_files", "pt_rotated_90.0.xyz")
        )

    def setUp(self, lpad=False):
        super().setUp(lpad=False)

    def tearDown(self):
        pass

    def test_rotate_torsion(self):
        atom_indexes = [6, 8, 9, 10]
        angle = 90.0
        ft = RotateTorsion(
            {"molecule": self.pt_mol, "atom_indexes": atom_indexes, "angle": angle}
        )
        rot_mol = ft.run_task({})
        test_mol = Molecule.from_dict(
            rot_mol.as_dict()["update_spec"]["prev_calc_molecule"]
        )
        np.testing.assert_equal(self.pt_rot_90_mol.species, test_mol.species)
        np.testing.assert_allclose(
            self.pt_rot_90_mol.cart_coords, test_mol.cart_coords, atol=0.0001
        )


class TestPerturbGeometry(AtomateTest):
    @classmethod
    def setUpClass(cls):

        cls.ts_init = Molecule.from_file(
            os.path.join(module_dir, "..", "..", "test_files", "ts_init.xyz")
        )
        cls.ts_perturbed = Molecule.from_file(
            os.path.join(module_dir, "..", "..", "test_files", "ts_perturbed.xyz")
        )
        cls.mode = QCOutput(
            os.path.join(module_dir, "..", "..", "test_files", "ts.out")
        ).data["frequency_mode_vectors"][0]

    def setUp(self, lpad=False):
        super().setUp(lpad=False)

    def tearDown(self):
        pass

    def test_perturb(self):
        ft = PerturbGeometry(
            {"molecule": self.ts_init, "mode": self.mode, "scale": 1.0}
        )
        pert_mol = ft.run_task({})
        test_mol = Molecule.from_dict(
            pert_mol.as_dict()["update_spec"]["prev_calc_molecule"]
        )
        np.testing.assert_equal(self.ts_perturbed.species, test_mol.species)
        np.testing.assert_allclose(
            self.ts_perturbed.cart_coords, test_mol.cart_coords, atol=0.0001
        )


if __name__ == "__main__":
    unittest.main()
