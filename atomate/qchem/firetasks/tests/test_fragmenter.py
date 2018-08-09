# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
from monty.serialization import loadfn#, dumpfn
try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import build_MoleculeGraph, MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN

from atomate.qchem.firetasks.fragmenter import FragmentMolecule, build_unique_molecules, build_new_FWs
from atomate.utils.testing import AtomateTest

import networkx as nx
import networkx.algorithms.isomorphism as iso


__author__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class TestFragmentMolecule(AtomateTest):
    @classmethod
    def setUpClass(cls):
        cls.pc = Molecule.from_file(
            os.path.join(module_dir, "..", "..", "test_files", "PC.xyz"))
        cls.pc_edges = [[5, 10], [5, 12], [5, 11], [5, 3], [3, 7], [3, 4],
                        [3, 0], [4, 8], [4, 9], [4, 1], [6, 1], [6, 0], [6, 2]]
        cls.pc_frag1 = Molecule.from_file(
            os.path.join(module_dir, "..", "..", "test_files", "PC_frag1.xyz"))
        cls.pc_frag1_edges = [[0, 2], [4, 2], [2, 1], [1, 3]]
        cls.tfsi = Molecule.from_file(os.path.join(module_dir, "..", "..", "test_files", "TFSI.xyz"))
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

    def test_build_unique_molecules(self):
        mol_graph = build_MoleculeGraph(self.pc, edges=self.pc_edges)
        unique_fragments = mol_graph.build_unique_fragments()
        unique_molecules = build_unique_molecules(unique_fragments,
                                                  self.pc.charge)
        self.assertEqual(len(unique_molecules), 295 * 3)
        # dumpfn(unique_molecules, os.path.join(module_dir,"pc_mols.json"))
        ref_mols = loadfn(os.path.join(module_dir, "pc_mols.json"))
        self.assertEqual(unique_molecules, ref_mols)

        mol_graph = build_MoleculeGraph(self.pc_frag1, edges=self.pc_frag1_edges)
        unique_fragments = mol_graph.build_unique_fragments()
        unique_molecules = build_unique_molecules(unique_fragments,
                                                  self.pc.charge)
        self.assertEqual(len(unique_molecules), 12 * 3)
        # dumpfn(unique_molecules, os.path.join(module_dir,"pc_frag1_mols.json"))
        ref_mols = loadfn(os.path.join(module_dir, "pc_frag1_mols.json"))
        self.assertEqual(unique_molecules, ref_mols)

    def test_build_new_FWs(self):
        mol_graph = build_MoleculeGraph(self.pc_frag1, edges=self.pc_frag1_edges)
        unique_fragments = mol_graph.build_unique_fragments()
        unique_molecules = build_unique_molecules(unique_fragments,
                                                  self.pc.charge)
        new_FWs = build_new_FWs(unique_molecules, [], 32, {})
        self.assertEqual(len(new_FWs), 36)
        doc = loadfn(os.path.join(module_dir, "doc.json"))
        new_FWs = build_new_FWs(unique_molecules, doc, 32, {})
        self.assertEqual(len(new_FWs), 29)

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

    def test_fragmenter(self):
        nm = iso.categorical_node_match("specie", "ERROR")
        mol_names = ["BF4-.xyz","DEC.xyz","DMC.xyz","EC.xyz","EMC.xyz","FEC.xyz","FSI-.xyz","PC.xyz","PF6-.xyz","TFSI-.xyz","VC.xyz"]
        num_frags = [5, 316, 43, 69, 194, 133, 35, 295, 7, 156, 37]
        all_fragments = []
        for ii,name in enumerate(mol_names):
            mol = Molecule.from_file(os.path.join(module_dir, "..", "..", "test_files", "top_11", name))
            mol_graph = build_MoleculeGraph(mol, strategy=OpenBabelNN,
                                            reorder=False, extend_structure=False)
            unique_fragments = mol_graph.build_unique_fragments()
            self.assertEqual(len(unique_fragments),num_frags[ii])
            for fragment in unique_fragments:
                all_fragments.append(fragment)
        self.assertEqual(len(all_fragments),1290)
        for fragment in all_fragments:
            if not [nx.is_isomorphic(fragment, f, node_match=nm) for f in unique_fragments].count(True) >= 1:
                unique_fragments.append(fragment)
        self.assertEqual(len(unique_fragments),834)

if __name__ == "__main__":
    unittest.main()
