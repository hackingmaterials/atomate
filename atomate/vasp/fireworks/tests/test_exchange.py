# coding: utf-8

import os
import unittest
import pandas as pd

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.transformations.standard_transformations import RotationTransformation

from atomate.vasp.fireworks.exchange import *

__author__ = "Nathan C. Frey"
__email__ = "ncfrey@lbl.gov"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
test_dir = os.path.join(module_dir, "..", "..", "test_files", "exchange_wf")


class TestExchangeFireworks(unittest.TestCase):
    def setUp(self):
        self.Mn3Al = pd.read_json(os.path.join(test_dir, "Mn3Al.json"))
        self.db_file = ""
        self.uuid = 1
        self.structures = [Structure.from_dict(s) for s in self.Mn3Al.structure]
        self.parent_structure = self.structures[0]
        self.energies = [
            e * len(self.parent_structure) for e in self.Mn3Al.energy_per_atom
        ]
        self.heisenberg_settings = {"cutoff": 3.0, "tol": 0.04, "avg": False}
        self.db_file = os.path.join(db_dir, "db.json")

    def test_EFWs(self):
        hmfw = HeisenbergModelFW(
            exchange_wf_uuid=self.uuid,
            parent_structure=self.parent_structure,
            parents=None,
            heisenberg_settings=self.heisenberg_settings,
        )
        num_ftasks = len(hmfw.tasks)
        self.assertEqual(num_ftasks, 5)

        vcfw = VampireCallerFW(
            exchange_wf_uuid=self.uuid,
            parent_structure=self.parent_structure,
            parents=None,
        )
        num_ftasks = len(vcfw.tasks)
        self.assertEqual(num_ftasks, 2)


if __name__ == "__main__":
    unittest.main()
