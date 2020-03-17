# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
import pandas as pd

from monty.os.path import which

from atomate.vasp.workflows.base.exchange import ExchangeWF
from atomate.vasp.firetasks.parse_outputs import (
    MagneticDeformationToDB,
    MagneticOrderingsToDB,
)
from atomate.utils.testing import AtomateTest, DB_DIR

from json import load
from pymatgen import Structure

__author__ = "Nathan C. Frey"
__email__ = "ncfrey@lbl.gov"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
test_dir = os.path.join(module_dir, "..", "..", "test_files", "exchange_wf")

enum_cmd = which("enum.x") or which("multienum.x")
VAMPEXE = which("vampire-serial")
enumlib_present = enum_cmd
vampire_present = VAMPEXE


class TestExchangeWF(AtomateTest):
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
        cls.cutoff = 3.0
        cls.tol = 0.04
        cls.db_file = os.path.join(db_dir, "db.json")

    @unittest.skipIf(not enumlib_present, "enumlib not present")
    def test_workflow(self):

        wf = ExchangeWF(
            magnetic_structures=self.structures,
            energies=self.energies,
            db_file=self.db_file,
        ).get_wf()
        self.assertEqual(wf.name, "Mn3Al - Exchange")


if __name__ == "__main__":
    unittest.main()
