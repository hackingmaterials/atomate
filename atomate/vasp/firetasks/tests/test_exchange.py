# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
import pandas as pd

from fireworks.utilities.fw_serializers import load_object

from atomate.vasp.firetasks.exchange_tasks import (
    HeisenbergModelMapping,
    HeisenbergConvergence,
    VampireMC,
)

from atomate.utils.testing import AtomateTest

from pymatgen import Structure
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.vasp import Incar, Poscar, Potcar, Kpoints
from pymatgen.io.vasp.sets import MPRelaxSet

__author__ = "Nathan C. Frey"
__email__ = "ncfrey@lbl.gov"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
test_dir = os.path.join(module_dir, "..", "..", "test_files", "exchange_wf")
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")

DEBUG_MODE = (
    False
)  # If true, retains the database and output dirs at the end of the test
VASP_CMD = (
    None
)  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...


class TestExchangeTasks(AtomateTest):
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

        new_fw_spec = {"_fw_env": {"db_file": os.path.join(db_dir, "db.json")}}

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_heisenberg_mm(self):
        d = dict(
            db_file=self.db_file,
            exchange_wf_uuid=self.uuid,
            parent_structure=self.parent_structure,
            cutoff=self.cutoff,
            tol=self.tol,
        )
        hmm = HeisenbergModelMapping(d)
        hmm.run_task({})

    def test_heisenberg_convergence(self):
        d = dict(
            db_file=self.db_file,
            exchange_wf_uuid=self.uuid,
            parent_structure=self.parent_structure,
        )
        hc = HeisenbergConvergence(d)
        hc.run_task({})

    def test_vampire_mc(self):
        d = dict(
            db_file=self.db_file,
            exchange_wf_uuid=self.uuid,
            parent_structure=self.parent_structure,
        )
        vmc = VampireMC(d)
        vmc.run_task({})


if __name__ == "__main__":
    unittest.main()
