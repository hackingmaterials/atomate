# coding: utf-8
# Copyright (c) Materials Virtual Lab.
# Distributed under the terms of the BSD License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import os
import unittest

from pymatgen.io.vasp import Outcar, Oszicar

from atomate.vasp.drones import VaspDrone


module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class VaspToDbTaskDroneTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.relax = os.path.join(module_dir, "..", "test_files", "Si_structure_optimization",
                                 "outputs")
        cls.relax2 = os.path.join(module_dir, "..", "test_files", "Si_structure_optimization_relax2",
                                 "outputs")
        cls.Al = os.path.join(module_dir, "..", "test_files", "Al")

    def test_assimilate(self):
        drone = VaspDrone()
        doc = drone.assimilate(self.relax)
        # Only the main changes from the vasprun as dict format and currently
        # used schema in pymatgen-db are tested for now.
        self.assertEqual(doc["composition_reduced"], {'Si': 1.0})
        self.assertEqual(doc["composition_unit_cell"], {'Si': 2.0})
        self.assertAlmostEqual(doc["output"]["energy"], -10.84671647)
        self.assertEqual(doc["formula_pretty"], 'Si')
        self.assertEqual(doc["formula_anonymous"], 'A')
        self.assertEqual(doc["calcs_reversed"][0]["output"]["energy"], doc["output"]["energy"])
        self.assertEqual(doc["input"]["parameters"]["ISMEAR"], -5)

    def test_runs_assimilate(self):
        drone = VaspDrone(runs=["relax1", "relax2"])
        doc = drone.assimilate(self.relax2)
        oszicar2 = Oszicar(os.path.join(self.relax2, "OSZICAR.relax2.gz"))
        outcar1 = Outcar(os.path.join(self.relax2, "OUTCAR.relax1.gz"))
        outcar2 = Outcar(os.path.join(self.relax2, "OUTCAR.relax2.gz"))
        outcar1 = outcar1.as_dict()
        outcar2 = outcar2.as_dict()
        run_stats1 = outcar1.pop("run_stats")
        run_stats2 = outcar2.pop("run_stats")
        self.assertEqual(len(doc["calcs_reversed"]), 2)
        self.assertEqual(doc["composition_reduced"], {'Si': 1.0})
        self.assertEqual(doc["composition_unit_cell"], {'Si': 2.0})
        self.assertAlmostEqual(doc["output"]["energy"], oszicar2.ionic_steps[-1]["E0"])
        self.assertEqual(doc["formula_pretty"], 'Si')
        self.assertEqual(doc["formula_anonymous"], 'A')
        self.assertEqual(list(doc["calcs_reversed"][0]["input"].keys()),
                         list(doc["calcs_reversed"][1]["input"].keys()))
        self.assertEqual(list(doc["calcs_reversed"][0]["output"].keys()),
                         list(doc["calcs_reversed"][1]["output"].keys()))
        self.assertEqual(doc["calcs_reversed"][0]["output"]["energy"], doc["output"]["energy"])
        self.assertEqual(doc["run_stats"][doc["calcs_reversed"][0]["task"]["name"]], run_stats2)
        self.assertEqual(doc["run_stats"][doc["calcs_reversed"][1]["task"]["name"]], run_stats1)
        self.assertEqual(doc["calcs_reversed"][0]["output"]["outcar"], outcar2)
        self.assertEqual(doc["calcs_reversed"][1]["output"]["outcar"], outcar1)

    def test_bandstructure(self):
        drone = VaspDrone()
        doc = drone.assimilate(self.Al)
        self.assertEqual(doc["composition_reduced"], {'Al': 1.0})
        self.assertEqual(doc["formula_pretty"], 'Al')
        self.assertEqual(doc["formula_anonymous"], 'A')
        for d in [doc["calcs_reversed"][0]["output"], doc["output"]]:
            self.assertIsNone(d["vbm"])
            self.assertIsNone(d["cbm"])
            self.assertEqual(d["bandgap"], 0.0)
            self.assertFalse(d["is_gap_direct"])
            self.assertTrue(d["is_metal"])

if __name__ == "__main__":
    unittest.main()
