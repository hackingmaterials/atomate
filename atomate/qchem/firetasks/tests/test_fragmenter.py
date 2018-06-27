# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
import shutil
from monty.serialization import loadfn, dumpfn
try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from atomate.qchem.firetasks.fragmenter import FragmentMolecule, edges_from_babel, build_unique_fragments, build_unique_molecules, build_new_FWs, build_MoleculeGraph, is_isomorphic
from atomate.utils.testing import AtomateTest
from pymatgen.core import Molecule
from pymatgen.analysis.graphs import MoleculeGraph

__author__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class TestFragmentMolecule(AtomateTest):
    @classmethod
    def setUpClass(cls):
        cls.pc = Molecule.from_file(os.path.join(module_dir, "..", "..", "test_files", "PC.xyz"))
        cls.pc_edges = [[5,10],[5,12],[5,11],[5,3],[3,7],[3,4],[3,0],[4,8],[4,9],[4,1],[6,1],[6,0],[6,2]]
        cls.pc_frag1 = Molecule.from_file(os.path.join(module_dir, "..", "..", "test_files", "PC_frag1.xyz"))
        cls.pc_frag1_edges = [[0,2],[4,2],[2,1],[1,3]]

    def setUp(self, lpad=False):
        super(TestFragmentMolecule, self).setUp(lpad=False)

    def tearDown(self):
        pass

    def test_edges_given_PC(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction") as FWAction_patch:
            ft = FragmentMolecule(molecule=self.pc, edges=self.pc_edges)
            ft.run_task({})
            self.assertEqual(len(FWAction_patch.call_args[1]["additions"]),295*3)

    def test_babel_PC(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction") as FWAction_patch:
            ft = FragmentMolecule(molecule=self.pc)
            ft.run_task({})
            self.assertEqual(len(FWAction_patch.call_args[1]["additions"]),295*3)

    def test_edges_given_PC_frag1(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction") as FWAction_patch:
            ft = FragmentMolecule(molecule=self.pc_frag1, edges=self.pc_frag1_edges)
            ft.run_task({})
            self.assertEqual(len(FWAction_patch.call_args[1]["additions"]),12*3)
        
    def test_babel_PC_frag1(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction") as FWAction_patch:
            ft = FragmentMolecule(molecule=self.pc_frag1)
            ft.run_task({})
            self.assertEqual(len(FWAction_patch.call_args[1]["additions"]),12*3)

    def test_edges_from_babel(self):
        self.assertEqual(edges_from_babel(self.pc),self.pc_edges)
        self.assertEqual(edges_from_babel(self.pc_frag1, self.pc_frag1_edges))

    def test_build_MoleculeGraph(self):
        mol_graph = build_MoleculeGraph(self.pc_frag1, self.pc_frag1_edges)
        # dumpfn(mol_graph.as_dict(), os.path.join(os.path.dirname(__file__),"pc_frag1_mg.json"))
        ref_mol_graph = loadfn(os.path.join(os.path.dirname(__file__),"pc_frag1_mg.json"))
        self.assertEqual(mol_graph, ref_mol_graph)
        self.assertEqual(mol_graph.graph.adj, ref_mol_graph.graph.adj)
        for node in mol_graph.graph:
            self.assertEqual(mol_graph.graph.node[node]["specie"],ref_mol_graph.graph.node[node]["specie"])
            for ii in range(3):
                self.assertEqual(mol_graph.graph.node[node]["coords"][ii],ref_mol_graph.graph.node[node]["coords"]["data"][ii])

        mol_graph = build_MoleculeGraph(self.pc, self.pc_edges)
        # dumpfn(mol_graph.as_dict(), os.path.join(os.path.dirname(__file__),"pc_mg.json"))
        ref_mol_graph = loadfn(os.path.join(os.path.dirname(__file__),"pc_mg.json"))
        self.assertEqual(mol_graph, ref_mol_graph)
        self.assertEqual(mol_graph.graph.adj, ref_mol_graph.graph.adj)
        for node in mol_graph.graph:
            self.assertEqual(mol_graph.graph.node[node]["specie"],ref_mol_graph.graph.node[node]["specie"])
            for ii in range(3):
                self.assertEqual(mol_graph.graph.node[node]["coords"][ii],ref_mol_graph.graph.node[node]["coords"]["data"][ii])

    def test_build_unique_fragments(self):
        mol_graph = build_MoleculeGraph(self.pc, self.pc_edges)
        unique_fragments = build_unique_fragments(mol_graph)
        self.assertEqual(len(unique_fragments),295)
        for ii in range(295):
            for jj in range(ii+1,295):
                self.assertEqual(is_isomorphic(unique_fragments[ii],unique_fragments[jj]),False)

    def test_build_unique_molecules(self):
        mol_graph = build_MoleculeGraph(self.pc, self.pc_edges)
        unique_fragments = build_unique_fragments(mol_graph)
        unique_molecules = build_unique_molecules(unique_fragments, self.pc.charge)
        self.assertEqual(len(unique_molecules), 295*3)
        # dumpfn(unique_molecules, os.path.join(os.path.dirname(__file__),"pc_mols.json"))
        ref_mols = loadfn(os.path.join(os.path.dirname(__file__),"pc_mols.json"))
        self.assertEqual(unique_molecules, ref_mols)

        mol_graph = build_MoleculeGraph(self.pc_frag1, self.pc_frag1_edges)
        unique_fragments = build_unique_fragments(mol_graph)
        unique_molecules = build_unique_molecules(unique_fragments, self.pc.charge)
        self.assertEqual(len(unique_molecules), 12*3)
        # dumpfn(unique_molecules, os.path.join(os.path.dirname(__file__),"pc_frag1_mols.json"))
        ref_mols = loadfn(os.path.join(os.path.dirname(__file__),"pc_frag1_mols.json"))
        self.assertEqual(unique_molecules, ref_mols)

    def test_build_new_FWs(self):
        mol_graph = build_MoleculeGraph(self.pc_frag1, self.pc_frag1_edges)
        unique_fragments = build_unique_fragments(mol_graph)
        unique_molecules = build_unique_molecules(unique_fragments, self.pc.charge)
        new_FWs = build_new_FWs(unique_molecules, [], 32, {})
        self.assertEqual(len(new_FWs), 36)
        doc = loadfn(os.path.join(os.path.dirname(__file__),"doc.json"))
        new_FWs = build_new_FWs(unique_molecules, doc, 32, {})
        self.assertEqual(len(new_FWs), 29)

if __name__ == "__main__":
    unittest.main()
