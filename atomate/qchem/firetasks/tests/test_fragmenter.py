# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
from monty.serialization import loadfn, dumpfn
try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import build_MoleculeGraph

from atomate.qchem.firetasks.fragmenter import FragmentMolecule
from atomate.qchem.firetasks.parse_outputs import QChemToDb
from atomate.utils.testing import AtomateTest


__author__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")


class TestFragmentMolecule(AtomateTest):
    @classmethod
    def setUpClass(cls):
        cls.pc = Molecule.from_file(
            os.path.join(module_dir, "..", "..", "test_files", "PC.xyz"))
        cls.pos_pc = Molecule.from_file(
            os.path.join(module_dir, "..", "..", "test_files", "PC.xyz"))
        cls.pos_pc.set_charge_and_spin(charge=1)
        cls.pc_edges = [[5, 10], [5, 12], [5, 11], [5, 3], [3, 7], [3, 4],
                        [3, 0], [4, 8], [4, 9], [4, 1], [6, 1], [6, 0], [6, 2]]
        cls.pc_frag1 = Molecule.from_file(
            os.path.join(module_dir, "..", "..", "test_files", "PC_frag1.xyz"))
        cls.pc_frag1_edges = [[0, 2], [4, 2], [2, 1], [1, 3]]
        cls.tfsi = Molecule.from_file(os.path.join(module_dir, "..", "..", "test_files", "TFSI.xyz"))
        cls.neg_tfsi = Molecule.from_file(os.path.join(module_dir, "..", "..", "test_files", "TFSI.xyz"))
        cls.neg_tfsi.set_charge_and_spin(charge=-1)
        cls.tfsi_edges = [14,1],[1,4],[1,5],[1,7],[7,11],[7,12],[7,13],[14,0],[0,2],[0,3],[0,6],[6,8],[6,9],[6,10]

    def setUp(self, lpad=False):
        super(TestFragmentMolecule, self).setUp(lpad=False)

    def tearDown(self):
        pass

    def test_edges_given_PC(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(molecule=self.pc, edges=self.pc_edges)
            ft.run_task({})
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 295 * 3)

    def test_babel_PC(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(molecule=self.pc)
            ft.run_task({})
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 295 * 3)

    def test_edges_given_PC_frag1(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(
                molecule=self.pc_frag1, edges=self.pc_frag1_edges)
            ft.run_task({})
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 12 * 3)

    def test_babel_PC_frag1(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(molecule=self.pc_frag1)
            ft.run_task({})
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 12 * 3)

    def test_edges_given_TFSI(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(molecule=self.tfsi, edges=self.tfsi_edges)
            ft.run_task({})
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 468)

    def test_babel_TFSI(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(molecule=self.tfsi)
            ft.run_task({})
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 468)

    def test_neg_TFSI(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(molecule=self.neg_tfsi)
            ft.run_task({})
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 624)

    def test_build_unique_relevant_molecules(self):
        FM = FragmentMolecule(molecule=self.pc, edges=self.pc_edges)
        FM.mol = FM.get("molecule")
        mol_graph = build_MoleculeGraph(self.pc, edges=self.pc_edges)
        FM.unique_fragments = mol_graph.build_unique_fragments()
        FM._build_unique_relevant_molecules()
        self.assertEqual(len(FM.unique_molecules), 295 * 3)
        # dumpfn(FM.unique_molecules, os.path.join(module_dir,"pc_mols.json"))
        ref_mols = loadfn(os.path.join(module_dir, "pc_mols.json"))
        self.assertEqual(FM.unique_molecules, ref_mols)

        FM = FragmentMolecule(molecule=self.pos_pc, edges=self.pc_edges)
        FM.mol = FM.get("molecule")
        mol_graph = build_MoleculeGraph(self.pos_pc, edges=self.pc_edges)
        FM.unique_fragments = mol_graph.build_unique_fragments()
        FM._build_unique_relevant_molecules()
        self.assertEqual(len(FM.unique_molecules), 295 * 4)
        # dumpfn(FM.unique_molecules, os.path.join(module_dir,"pos_pc_mols.json"))
        ref_mols = loadfn(os.path.join(module_dir, "pos_pc_mols.json"))
        self.assertEqual(FM.unique_molecules, ref_mols)

        FM = FragmentMolecule(molecule=self.pc_frag1, edges=self.pc_frag1_edges)
        FM.mol = FM.get("molecule")
        mol_graph = build_MoleculeGraph(self.pc_frag1, edges=self.pc_frag1_edges)
        FM.unique_fragments = mol_graph.build_unique_fragments()
        FM._build_unique_relevant_molecules()
        self.assertEqual(len(FM.unique_molecules), 12 * 3)
        # dumpfn(FM.unique_molecules, os.path.join(module_dir,"pc_frag1_mols.json"))
        ref_mols = loadfn(os.path.join(module_dir, "pc_frag1_mols.json"))
        self.assertEqual(FM.unique_molecules, ref_mols)

    def test_build_new_FWs(self):
        FM = FragmentMolecule(molecule=self.pc_frag1, edges=self.pc_frag1_edges)
        FM.mol = FM.get("molecule")
        mol_graph = build_MoleculeGraph(self.pc_frag1, edges=self.pc_frag1_edges)
        FM.unique_fragments = mol_graph.build_unique_fragments()
        FM._build_unique_relevant_molecules()
        FM.all_relevant_docs = list()
        new_FWs = FM._build_new_FWs()
        self.assertEqual(len(new_FWs), 36)

    def test_in_database_through_build_new_FWs(self):
        FM = FragmentMolecule(molecule=self.pc_frag1, edges=self.pc_frag1_edges)
        FM.mol = FM.get("molecule")
        mol_graph = build_MoleculeGraph(self.pc_frag1, edges=self.pc_frag1_edges)
        FM.unique_fragments = mol_graph.build_unique_fragments()
        FM._build_unique_relevant_molecules()
        FM.all_relevant_docs = loadfn(os.path.join(module_dir, "doc.json"))
        new_FWs = FM._build_new_FWs()
        self.assertEqual(len(new_FWs), 29)

    def test_in_database_with_actual_database(self):
        db_file=os.path.join(db_dir, "db.json")
        parse_firetask = QChemToDb(calc_dir=os.path.join(module_dir, "..", "..", "test_files", "2620_complete"), db_file=db_file)
        parse_firetask.run_task({})
        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(molecule=self.neg_tfsi, db_file=db_file)
            ft.run_task({})
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 623)


if __name__ == "__main__":
    unittest.main()
