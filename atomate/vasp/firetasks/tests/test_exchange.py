# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
import pandas as pd

from monty.os.path import which

from fireworks.utilities.fw_serializers import load_object

from atomate.vasp.firetasks.exchange import (
    HeisenbergModelMapping,
    HeisenbergModelToDb,
    HeisenbergConvergence,
    VampireMC,
    VampireToDb,
)

from atomate.utils.testing import AtomateTest

from pymatgen.core import Structure
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.vasp import Incar, Poscar, Potcar, Kpoints
from pymatgen.io.vasp.sets import MPRelaxSet

__author__ = "Nathan C. Frey"
__email__ = "ncfrey@lbl.gov"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
test_dir = os.path.join(module_dir, "..", "..", "test_files", "exchange_wf")
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")

VAMPEXE = which("vampire-serial")
vampire_present = VAMPEXE


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
        cls.heisenberg_settings = {"cutoff": 3.0, "tol": 0.04}
        cls.mc_settings = {
            "mc_box_size": 3,
            "equil_timesteps": 10,
            "mc_timesteps": 10,
            "avg": True,
        }

        cls.db_file = os.path.join(db_dir, "db.json")

        new_fw_spec = {"_fw_env": {"db_file": os.path.join(db_dir, "db.json")}}

    @unittest.skipIf(not vampire_present, "vampire not present")
    def test_heisenberg_mm(self):
        hmm = HeisenbergModelMapping(
            structures=self.structures,
            energies=self.energies,
            heisenberg_settings=self.heisenberg_settings,
        )
        hmm.run_task({})

        hmtdb = HeisenbergModelToDb(db_file=self.db_file, wf_uuid=self.uuid)
        hmtdb.run_task({})

        vmc = VampireMC(
            db_file=self.db_file, wf_uuid=self.uuid, mc_settings=self.mc_settings
        )
        vmc.run_task({})

        vtdb = VampireToDb(db_file=self.db_file, wf_uuid=self.uuid)
        vtdb.run_task({})


if __name__ == "__main__":
    unittest.main()
