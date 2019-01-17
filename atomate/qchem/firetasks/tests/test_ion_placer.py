# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
import json
try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from pymatgen.core.structure import Molecule
from pymatgen.io.qchem.outputs import QCOutput

from atomate.qchem.firetasks.ion_placer import PlaceIon
from atomate.qchem.firetasks.parse_outputs import QChemToDb
from atomate.qchem.database import QChemCalcDb
from atomate.utils.testing import AtomateTest


__author__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")


class TestPlaceIon(AtomateTest):
    @classmethod
    def setUpClass(cls):
        cls.pc_neut1_data = QCOutput(os.path.join(module_dir, "..", "..", "test_files", "ion_placer_files","PC01.qout")).data
        cls.ec_neg2_data = QCOutput(os.path.join(module_dir, "..", "..", "test_files", "ion_placer_files","EC-12.qout")).data

    def setUp(self, lpad=False):
        super(TestPlaceIon, self).setUp(lpad=False)

    def tearDown(self):
        pass

    def test_EC_neg_doublet(self):
        with patch("atomate.qchem.firetasks.ion_placer.FWAction") as FWAction_patch:
            ft = PlaceIon(molecule=self.ec_neg2_data["initial_molecule"], mulliken=self.ec_neg2_data["Mulliken"])
            ft.run_task({})
            self.assertEqual(len(ft.ion_positions),len(FWAction_patch.call_args[1]["detours"]))

        mol = self.ec_neg2_data["initial_molecule"]
        mol.add_site_property("charge",self.ec_neg2_data["Mulliken"][0][::,0])
        mol.add_site_property("spin",self.ec_neg2_data["Mulliken"][0][::,1])

        with patch("atomate.qchem.firetasks.ion_placer.FWAction") as FWAction_patch:
            ft = PlaceIon(molecule=mol)
            ft.run_task({})
            self.assertEqual(len(ft.ion_positions),len(FWAction_patch.call_args[1]["detours"]))

        with patch("atomate.qchem.firetasks.ion_placer.FWAction") as FWAction_patch:
            ft = PlaceIon(molecule=mol, charges=[-1,0,1])
            ft.run_task({})
            self.assertEqual(len(ft.ion_positions)*5,len(FWAction_patch.call_args[1]["detours"]))


    def test_PC_neut_singlet(self):
        with patch("atomate.qchem.firetasks.ion_placer.FWAction") as FWAction_patch:
            ft = PlaceIon(molecule=self.pc_neut1_data["initial_molecule"], mulliken=self.pc_neut1_data["Mulliken"])
            ft.run_task({})
            self.assertEqual(len(ft.ion_positions),len(FWAction_patch.call_args[1]["detours"]))

        mol = self.pc_neut1_data["initial_molecule"]
        mol.add_site_property("charge",self.pc_neut1_data["Mulliken"][0])

        with patch("atomate.qchem.firetasks.ion_placer.FWAction") as FWAction_patch:
            ft = PlaceIon(molecule=mol)
            ft.run_task({})
            self.assertEqual(len(ft.ion_positions),len(FWAction_patch.call_args[1]["detours"]))

        with patch("atomate.qchem.firetasks.ion_placer.FWAction") as FWAction_patch:
            ft = PlaceIon(molecule=mol, charges=[-1,0,1])
            ft.run_task({})
            self.assertEqual(len(ft.ion_positions)*5,len(FWAction_patch.call_args[1]["detours"]))


if __name__ == "__main__":
    unittest.main()
