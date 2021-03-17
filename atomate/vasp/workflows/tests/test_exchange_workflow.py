# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
import pandas as pd

from monty.os.path import which

from fireworks import FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.vasp.workflows.base.exchange import ExchangeWF
from atomate.vasp.firetasks.parse_outputs import (
    MagneticDeformationToDb,
    MagneticOrderingsToDb,
)
from atomate.vasp.database import VaspCalcDb
from atomate.utils.testing import AtomateTest, DB_DIR
from atomate.utils.utils import get_a_unique_id

from json import load
from pymatgen.core import Structure

__author__ = "Nathan C. Frey"
__email__ = "ncfrey@lbl.gov"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
test_dir = os.path.join(module_dir, "..", "..", "test_files", "exchange_wf")

VAMPEXE = which("vampire-serial")
vampire_present = VAMPEXE


class TestExchangeWF(AtomateTest):
    @classmethod
    def setUpClass(cls):
        cls.Mn3Al = pd.read_json(os.path.join(test_dir, "Mn3Al.json"))
        cls.wfid = "wfid_" + get_a_unique_id()
        cls.db_file = os.path.join(db_dir, "db.json")

        cls.structures = [Structure.from_dict(s) for s in cls.Mn3Al.structure]
        cls.parent_structure = cls.structures[0]
        cls.energies = [
            e * len(cls.parent_structure) for e in cls.Mn3Al.energy_per_atom
        ]
        cls.cutoff = 3.0
        cls.tol = 0.04
        cls.heisenberg_settings = {"cutoff": 3.0, "tol": 0.04}
        cls.mc_settings = {
            "mc_box_size": 3,
            "equil_timesteps": 10,
            "mc_timesteps": 10,
            "avg": False,
        }

    @unittest.skipIf(not vampire_present, "vampire not present")
    def test_workflow(self):
        wf = self._get_simulated_wflow()

        self.assertEqual(wf.name, "Mn3Al - Exchange")

        # Heisenberg + VampireMC = 2
        self.assertEqual(len(wf.fws), 2)

        fw_ids = self.lp.add_wf(wf)
        rapidfire(self.lp, fworker=FWorker(env={"db_file": self.db_file}))

        # Exchange collection
        exchange_collection = self.get_task_database().exchange

        # Check Heisenberg model mapping
        d = exchange_collection.find_one({"task_name": "heisenberg model"})
        self._check_run(d, "heisenberg model")

        # Check Vampire MC
        d = exchange_collection.find_one({"task_name": "vampire caller"})
        self._check_run(d, "vampire caller")

        # Check for "COMPLETED" status
        fw_id = list(fw_ids.values())[0]
        wf = self.lp.get_wf_by_fw_id(fw_id)
        is_completed = [s == "COMPLETED" for s in wf.fw_states.values()]
        
        self.assertTrue(all(is_completed))

    def _check_run(self, d, task_name):
        if task_name == "heisenberg model":
            self.assertEqual(d["nn_cutoff"], 3.0)
            self.assertAlmostEqual(d["heisenberg_model"]["javg"], 420.4, 2)

        if task_name == "vampire caller":
            self.assertAlmostEqual(d["vampire_output"]["critical_temp"], 25, 2)

    def _get_simulated_wflow(self):
        c = {}
        c["heisenberg_settings"] = self.heisenberg_settings
        c["mc_settings"] = self.mc_settings
        c["DB_FILE"] = self.db_file

        wf = ExchangeWF(
            magnetic_structures=self.structures,
            energies=self.energies,
            db_file=self.db_file,
        ).get_wf(c=c)

        return wf


if __name__ == "__main__":
    unittest.main()
