import os
import unittest

import pandas as pd
from pymatgen.core.structure import Structure

from atomate.utils.testing import AtomateTest
from atomate.vasp.fireworks.exchange import HeisenbergModelFW, VampireCallerFW

__author__ = "Nathan C. Frey"
__email__ = "ncfrey@lbl.gov"

module_dir = os.path.dirname(os.path.abspath(__file__))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
test_dir = os.path.join(module_dir, "..", "..", "test_files", "exchange_wf")


class TestExchangeFireworks(AtomateTest):
    @classmethod
    def setUpClass(cls):
        cls.Mn3Al = pd.read_json(os.path.join(test_dir, "Mn3Al.json"))
        cls.db_file = ""
        cls.uuid = 1
        cls.structures = [Structure.from_dict(s) for s in cls.Mn3Al.structure]
        cls.parent_structure = cls.structures[0]
        cls.energies = [
            e * len(cls.parent_structure) for e in cls.Mn3Al.energy_per_atom
        ]
        cls.heisenberg_settings = {"cutoff": 3.0, "tol": 0.04}
        cls.db_file = os.path.join(db_dir, "db.json")

    def test_EFWs(self):
        hmfw = HeisenbergModelFW(
            wf_uuid=self.uuid,
            parent_structure=self.parent_structure,
            parents=None,
            heisenberg_settings=self.heisenberg_settings,
            structures=self.structures,
            energies=self.energies,
        )
        num_ftasks = len(hmfw.tasks)
        self.assertEqual(num_ftasks, 2)

        vcfw = VampireCallerFW(
            wf_uuid=self.uuid, parent_structure=self.parent_structure, parents=[hmfw]
        )
        num_ftasks = len(vcfw.tasks)
        self.assertEqual(num_ftasks, 2)


if __name__ == "__main__":
    unittest.main()
