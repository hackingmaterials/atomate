# coding: utf-8
# Copyright (c) Materials Virtual Lab.
# Distributed under the terms of the BSD License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import os
import unittest

from pymatgen.io.vasp import Outcar, Oszicar

from atomate.vasp.drones import VaspDrone

import numpy as np


module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class VaspToDbTaskDroneTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.relax = os.path.join(module_dir, "..", "test_files", "Si_structure_optimization",
                                 "outputs")
        cls.relax2 = os.path.join(module_dir, "..", "test_files", "Si_structure_optimization_relax2",
                                 "outputs")
        cls.Al = os.path.join(module_dir, "..", "test_files", "Al")
        cls.Si_static = os.path.join(module_dir, "..", "test_files", "Si_static", "outputs")

    def test_assimilate(self):
        drone = VaspDrone()
        doc = drone.assimilate(self.relax)
        # Only the main changes from the vasprun as dict format and currently
        # used schema in pymatgen-db are tested for now.
        self.assertEqual(doc["composition_reduced"], {'Si': 1.0})
        self.assertEqual(doc["composition_unit_cell"], {'Si': 2.0})
        self.assertAlmostEqual(doc["output"]["energy"], -10.84671647)
        self.assertTrue(np.allclose(doc["output"]["forces"], [[0, 0, 0], [0, 0, 0]]))
        self.assertAlmostEqual(doc['output']['stress'][0][0], -0.08173155)
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

        doc = drone.assimilate(self.relax2)
        self.assertEqual(doc["composition_reduced"], {'Si': 1.0})
        self.assertEqual(doc["formula_pretty"], 'Si')
        self.assertEqual(doc["formula_anonymous"], 'A')
        for d in [doc["calcs_reversed"][0]["output"], doc["output"]]:
            self.assertAlmostEqual(d["vbm"],5.6147)
            self.assertAlmostEqual(d["cbm"],6.2652)
            self.assertAlmostEqual(d["bandgap"], 0.6505)
            self.assertFalse(d["is_gap_direct"])
            self.assertFalse(d["is_metal"])
            self.assertNotIn("transition",d)
            self.assertAlmostEqual(d["direct_gap"],2.5561)
            self.assertNotIn("bandstructure",doc["calcs_reversed"][0])


        doc = drone.assimilate(self.Si_static)
        self.assertEqual(doc["composition_reduced"], {'Si': 1.0})
        self.assertEqual(doc["formula_pretty"], 'Si')
        self.assertEqual(doc["formula_anonymous"], 'A')
        for d in [doc["calcs_reversed"][0]["output"], doc["output"]]:
            self.assertAlmostEqual(d["vbm"],5.6138)
            self.assertAlmostEqual(d["cbm"],6.2644)
            self.assertAlmostEqual(d["bandgap"],  0.6506)
            self.assertFalse(d["is_gap_direct"])
            self.assertFalse(d["is_metal"])
            self.assertNotIn("transition",d)
            self.assertAlmostEqual(d["direct_gap"],2.5563)
            self.assertIn("bandstructure", doc["calcs_reversed"][0])


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
            self.assertEqual(doc["calcs_reversed"][0]["bandstructure"]["@class"],"BandStructureSymmLine")

    def test_detect_output_file_paths(self):
        drone = VaspDrone()
        doc = drone.assimilate(self.Si_static)

        self.assertDictEqual({
            'chgcar': 'CHGCAR.gz',
            'locpot': 'LOCPOT.gz',
            'aeccar0': 'AECCAR0.gz',
            'aeccar1': 'AECCAR1.gz',
            'aeccar2': 'AECCAR2.gz',
            'procar': 'PROCAR.gz',
            'wavecar': 'WAVECAR.gz'
        }, doc['calcs_reversed'][0]['output_file_paths'])

        doc = drone.assimilate(self.relax2)
        self.assertDictEqual({'chgcar': 'CHGCAR.relax1.gz', 'procar': 'PROCAR.relax1.gz',
                              'wavecar': 'WAVECAR.relax1.gz'},
                             doc['calcs_reversed'][1]['output_file_paths'])

    def test_parse_locpot(self):
        drone = VaspDrone(parse_locpot=True)
        doc = drone.assimilate(self.Si_static)

        self.assertTrue(drone.parse_locpot)
        self.assertTrue('locpot' in doc['calcs_reversed'][0]['output'])
        self.assertTrue(0 in doc['calcs_reversed'][0]['output']['locpot'])
        self.assertTrue(1 in doc['calcs_reversed'][0]['output']['locpot'])
        self.assertTrue(2 in doc['calcs_reversed'][0]['output']['locpot'])

        self.assertAlmostEqual(np.sum(doc['calcs_reversed'][0]['output']['locpot'][0]),0)
        self.assertAlmostEqual(np.sum(doc['calcs_reversed'][0]['output']['locpot'][1]),0)
        self.assertAlmostEqual(np.sum(doc['calcs_reversed'][0]['output']['locpot'][2]),0)

    def test_parse_chrgcar(self):
        drone = VaspDrone(parse_chgcar=True, parse_aeccar=True)
        doc = drone.assimilate(self.Si_static)
        cc = doc['calcs_reversed'][0]['chgcar']
        self.assertAlmostEqual(cc.data['total'].sum()/cc.ngridpts, 8.0, 4)
        cc = doc['calcs_reversed'][0]['aeccar0']
        self.assertAlmostEqual(cc.data['total'].sum()/cc.ngridpts, 23.253588293583313, 4)
        cc = doc['calcs_reversed'][0]['aeccar2']
        self.assertAlmostEqual(cc.data['total'].sum()/cc.ngridpts, 8.01314480789829, 4)


if __name__ == "__main__":
    unittest.main()
