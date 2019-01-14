# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
import json
from monty.serialization import loadfn#, dumpfn
try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.io.qchem.outputs import QCOutput

from atomate.qchem.firetasks.fragmenter import FragmentMolecule
from atomate.qchem.firetasks.parse_outputs import QChemToDb
from atomate.qchem.database import QChemCalcDb
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
        cls.neg_ec = Molecule.from_file(
            os.path.join(module_dir, "..", "..", "test_files", "top_11", "EC.xyz"))
        cls.neg_ec.set_charge_and_spin(charge=-1)

    def setUp(self, lpad=False):
        super(TestFragmentMolecule, self).setUp(lpad=False)

    def tearDown(self):
        pass

    def test_edges_given_PC(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(molecule=self.pc, edges=self.pc_edges, depth=0, open_rings=False, do_triplets=False)
            ft.run_task({})
            self.assertEqual(ft.check_db,False)
            self.assertEqual(ft.charges,[-1,0,1])
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 295 * 3)

    def test_edges_given_PC_frag1(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(
                molecule=self.pc_frag1, edges=self.pc_frag1_edges, depth=0, do_triplets=False)
            ft.run_task({})
            self.assertEqual(ft.check_db,False)
            self.assertEqual(ft.charges,[-1,0,1])
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 12 * 3)

    def test_babel_PC_frag1(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(molecule=self.pc_frag1, depth=0, do_triplets=False)
            ft.run_task({})
            self.assertEqual(ft.check_db,False)
            self.assertEqual(ft.charges,[-1,0,1])
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 12 * 3)

    def test_edges_given_TFSI(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(molecule=self.tfsi, edges=self.tfsi_edges, depth=0, do_triplets=False)
            ft.run_task({})
            self.assertEqual(ft.check_db,False)
            self.assertEqual(ft.charges,[-1,0,1])
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 468)

    def test_babel_TFSI(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(molecule=self.tfsi, depth=0, do_triplets=False)
            ft.run_task({})
            self.assertEqual(ft.check_db,False)
            self.assertEqual(ft.charges,[-1,0,1])
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 468)

    def test_neg_TFSI_with_additional_charges(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(molecule=self.neg_tfsi, depth=0, additional_charges=[-2,1], do_triplets=False)
            ft.run_task({})
            self.assertEqual(ft.check_db,False)
            self.assertEqual(ft.charges,[-1,0,-2,1])
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 624)

        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(molecule=self.neg_tfsi, depth=0, additional_charges=[0,1], do_triplets=False)
            ft.run_task({})
            self.assertEqual(ft.check_db,False)
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 468)

    def test_neg_TFSI(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(molecule=self.neg_tfsi, depth=0, do_triplets=False)
            ft.run_task({})
            self.assertEqual(ft.check_db,False)
            self.assertEqual(ft.charges,[-1,0])
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 312)

    def test_build_unique_relevant_molecules(self):
        ft = FragmentMolecule(molecule=self.pc, edges=self.pc_edges, depth=0, opt_steps=1000)
        ft.mol = ft.get("molecule")
        ft.depth = ft.get("depth")
        ft.charges = [-1, 0, 1]
        ft.do_triplets = False
        edges = {(e[0], e[1]): None for e in self.pc_edges}
        mol_graph = MoleculeGraph.with_edges(self.pc, edges)
        ft.unique_fragments = mol_graph.build_unique_fragments()
        ft._build_unique_relevant_molecules()
        self.assertEqual(len(ft.unique_molecules), 295 * 3)
        #dumpfn(ft.unique_molecules, os.path.join(module_dir,"pc_mols.json"))
        ref_mols = loadfn(os.path.join(module_dir, "pc_mols.json"))
        self.assertEqual(ft.unique_molecules, ref_mols)

        ft = FragmentMolecule(molecule=self.pos_pc, edges=self.pc_edges, depth=0, opt_steps=1000)
        ft.mol = ft.get("molecule")
        ft.depth = ft.get("depth")
        ft.charges = [-1, 0, 1, 2]
        ft.do_triplets = False
        mol_graph = MoleculeGraph.with_edges(self.pos_pc, edges)
        ft.unique_fragments = mol_graph.build_unique_fragments()
        ft._build_unique_relevant_molecules()
        self.assertEqual(len(ft.unique_molecules), 295 * 4)
        #dumpfn(ft.unique_molecules, os.path.join(module_dir,"pos_pc_mols.json"))
        ref_mols = loadfn(os.path.join(module_dir, "pos_pc_mols.json"))
        self.assertEqual(ft.unique_molecules, ref_mols)

        ft = FragmentMolecule(molecule=self.pc_frag1, edges=self.pc_frag1_edges, depth=0)
        ft.mol = ft.get("molecule")
        ft.depth = ft.get("depth")
        ft.charges = [-1, 0, 1]
        ft.do_triplets = False
        pc_frag1_edges = {(e[0], e[1]): None for e in self.pc_frag1_edges}
        mol_graph = MoleculeGraph.with_edges(self.pc_frag1, pc_frag1_edges)
        ft.unique_fragments = mol_graph.build_unique_fragments()
        ft._build_unique_relevant_molecules()
        self.assertEqual(len(ft.unique_molecules), 12 * 3)
        #dumpfn(ft.unique_molecules, os.path.join(module_dir,"pc_frag1_mols.json"))
        ref_mols = loadfn(os.path.join(module_dir, "pc_frag1_mols.json"))
        self.assertEqual(ft.unique_molecules, ref_mols)

    def test_build_unique_relevant_molecules_with_triplets(self):
        ft = FragmentMolecule(molecule=self.pc, edges=self.pc_edges, depth=0, opt_steps=1000)
        ft.mol = ft.get("molecule")
        ft.depth = ft.get("depth")
        ft.charges = [-1, 0, 1]
        ft.do_triplets = True
        edges = {(e[0], e[1]): None for e in self.pc_edges}
        mol_graph = MoleculeGraph.with_edges(self.pc, edges)
        ft.unique_fragments = mol_graph.build_unique_fragments()
        ft._build_unique_relevant_molecules()
        self.assertEqual(len(ft.unique_molecules), 1323)
        #dumpfn(ft.unique_molecules, os.path.join(module_dir,"pc_mols_with_trips.json"))
        ref_mols = loadfn(os.path.join(module_dir, "pc_mols_with_trips.json"))
        self.assertEqual(ft.unique_molecules, ref_mols)

        ft = FragmentMolecule(molecule=self.pos_pc, edges=self.pc_edges, depth=0, opt_steps=1000)
        ft.mol = ft.get("molecule")
        ft.depth = ft.get("depth")
        ft.charges = [-1, 0, 1, 2]
        ft.do_triplets = True
        mol_graph = MoleculeGraph.with_edges(self.pos_pc, edges)
        ft.unique_fragments = mol_graph.build_unique_fragments()
        ft._build_unique_relevant_molecules()
        self.assertEqual(len(ft.unique_molecules), 1770)
        #dumpfn(ft.unique_molecules, os.path.join(module_dir,"pos_pc_mols_with_trips.json"))
        ref_mols = loadfn(os.path.join(module_dir, "pos_pc_mols_with_trips.json"))
        self.assertEqual(ft.unique_molecules, ref_mols)

        ft = FragmentMolecule(molecule=self.pc_frag1, edges=self.pc_frag1_edges, depth=0)
        ft.mol = ft.get("molecule")
        ft.depth = ft.get("depth")
        ft.charges = [-1, 0, 1]
        ft.do_triplets = True
        pc_frag1_edges = {(e[0], e[1]): None for e in self.pc_frag1_edges}
        mol_graph = MoleculeGraph.with_edges(self.pc_frag1, pc_frag1_edges)
        ft.unique_fragments = mol_graph.build_unique_fragments()
        ft._build_unique_relevant_molecules()
        self.assertEqual(len(ft.unique_molecules), 54)
        #dumpfn(ft.unique_molecules, os.path.join(module_dir,"pc_frag1_mols_with_trips.json"))
        ref_mols = loadfn(os.path.join(module_dir, "pc_frag1_mols_with_trips.json"))
        self.assertEqual(ft.unique_molecules, ref_mols)

    def test_build_new_FWs(self):
        ft = FragmentMolecule(molecule=self.pc_frag1, edges=self.pc_frag1_edges, depth=0)
        ft.mol = ft.get("molecule")
        ft.depth = ft.get("depth")
        ft.charges = [-1, 0, 1]
        ft.do_triplets = False
        ft.qchem_input_params = {}
        pc_frag1_edges = {(e[0], e[1]): None for e in self.pc_frag1_edges}
        mol_graph = MoleculeGraph.with_edges(self.pc_frag1, pc_frag1_edges)
        ft.unique_fragments = mol_graph.build_unique_fragments()
        ft._build_unique_relevant_molecules()
        ft.all_relevant_docs = list()
        new_FWs = ft._build_new_FWs()
        self.assertEqual(len(new_FWs), 36)

    def test_in_database_through_build_new_FWs(self):
        ft = FragmentMolecule(molecule=self.pc_frag1, edges=self.pc_frag1_edges, depth=0)
        ft.mol = ft.get("molecule")
        ft.depth = ft.get("depth")
        ft.charges = [-1, 0, 1]
        ft.do_triplets = False
        ft.qchem_input_params = {}
        pc_frag1_edges = {(e[0], e[1]): None for e in self.pc_frag1_edges}
        mol_graph = MoleculeGraph.with_edges(self.pc_frag1, pc_frag1_edges)
        ft.unique_fragments = mol_graph.build_unique_fragments()
        ft._build_unique_relevant_molecules()
        docs = loadfn(os.path.join(module_dir, "doc.json"))
        for doc in docs:
            doc["input"]["initial_molecule"] = doc["input"]["initial_molecule"].as_dict()
        ft.all_relevant_docs = docs
        new_FWs = ft._build_new_FWs()
        self.assertEqual(len(new_FWs), 29)

    def test_in_database_and_EC_neg_frag(self):
        db_file = os.path.join(db_dir, "db.json")
        mmdb = QChemCalcDb.from_db_file(db_file, admin=True)
        with open(os.path.join(module_dir, "..", "..", "test_files","sb40.json")) as f:
            tmp = json.load(f)
            for entry in tmp:
                mmdb.insert(entry)
        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(molecule=self.neg_ec, depth=1, qchem_input_params={"pcm_dielectric": 40.0}, check_db=True, db_file=db_file)
            ft.run_task({})
            self.assertEqual(ft.check_db,True)
            frags = ft.unique_fragments
            self.assertEqual(len(frags), 7)
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 0)
        mmdb.reset()

    def test_in_database_with_actual_database(self):
        db_file=os.path.join(db_dir, "db.json")
        dir2620=os.path.join(module_dir, "..", "..", "test_files", "2620_complete")
        mol2620=QCOutput(os.path.join(dir2620,"mol.qout.opt_0")).data["initial_molecule"]
        parse_firetask = QChemToDb(calc_dir=dir2620, db_file=db_file)
        parse_firetask.run_task({})
        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(molecule=self.neg_tfsi, db_file=db_file, depth=0, additional_charges=[-2,1], do_triplets=False)
            ft.run_task({})
            self.assertEqual(ft.check_db,True)
            self.assertEqual(ft.charges,[-1,0,-2,1])
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 623)
            self.assertEqual(ft._in_database(mol2620),True)
            mol2620.set_charge_and_spin(charge=0)
            self.assertEqual(ft._in_database(mol2620),False)

    def test_babel_PC_with_RO_depth_0_vs_depth_10(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(molecule=self.pc, depth=0, open_rings=True, do_triplets=False, opt_steps=1000)
            ft.run_task({})
            self.assertEqual(ft.check_db,False)
            depth0frags = ft.unique_fragments
            self.assertEqual(len(depth0frags), 509)
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 1527)

        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(molecule=self.pc, depth=10, open_rings=True, do_triplets=False, opt_steps=1000)
            ft.run_task({})
            self.assertEqual(ft.check_db,False)
            depth10frags = ft.unique_fragments
            self.assertEqual(len(depth10frags), 509)
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 1513)

        for fragment10 in depth10frags:
            found = False
            for fragment0 in depth0frags:
                if fragment0.isomorphic_to(fragment10):
                    found = True
                    break
            self.assertEqual(found, True)

    def test_babel_PC_depth_0_vs_depth_10(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(molecule=self.pc, depth=0, open_rings=False, do_triplets=False)
            ft.run_task({})
            self.assertEqual(ft.check_db,False)
            depth0frags = ft.unique_fragments
            self.assertEqual(len(depth0frags), 295)
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 295*3)

        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(molecule=self.pc, depth=10, open_rings=False, do_triplets=False)
            ft.run_task({})
            self.assertEqual(ft.check_db,False)
            depth10frags = ft.unique_fragments
            self.assertEqual(len(depth10frags), 63)
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 63*3)

        for fragment10 in depth10frags:
            found = False
            for fragment0 in depth0frags:
                if fragment0.isomorphic_to(fragment10):
                    found = True
            self.assertEqual(found, True)

    def test_EC_neg(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction"
                   ) as FWAction_patch:
            ft = FragmentMolecule(molecule=self.neg_ec, depth=1, check_db=False)
            ft.run_task({})
            self.assertEqual(ft.check_db,False)
            frags = ft.unique_fragments
            self.assertEqual(len(frags), 7)
            self.assertEqual(
                len(FWAction_patch.call_args[1]["additions"]), 15)


if __name__ == "__main__":
    unittest.main()
