# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
import shutil
try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from atomate.qchem.firetasks.fragmenter import FragmentMolecule
from atomate.utils.testing import AtomateTest
from pymatgen.core import Molecule
from pymatgen.analysis.graphs import MoleculeGraph

__author__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class TestFragmentMolecule(AtomateTest):

    def setUp(self, lpad=False):
        super(TestFragmentMolecule, self).setUp(lpad=False)

    def tearDown(self):
        pass

    def test_edges_given_PC(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction") as FWAction_patch:
            pc = Molecule.from_file(os.path.join(module_dir, "..", "..", "test_files", "PC.xyz"))
            edges = [[5,10],[5,12],[5,11],[5,3],[3,7],[3,4],[3,0],[4,8],[4,9],[4,1],[6,1],[6,0],[6,2]]
            ft = FragmentMolecule(molecule=pc, edges=edges)
            ft.run_task({})
            self.assertEqual(len(FWAction_patch.call_args[1]["additions"]),295)

    def test_babel_PC(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction") as FWAction_patch:
            pc = Molecule.from_file(os.path.join(module_dir, "..", "..", "test_files", "PC.xyz"))
            ft = FragmentMolecule(molecule=pc)
            ft.run_task({})
            self.assertEqual(len(FWAction_patch.call_args[1]["additions"]),295)

    def test_edges_given_PC_frag1(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction") as FWAction_patch:
            pc_frag1 = Molecule.from_file(os.path.join(module_dir, "..", "..", "test_files", "PC_frag1.xyz"))
            edges = [[0,2],[4,2],[2,1],[1,3]]
            ft = FragmentMolecule(molecule=pc_frag1, edges=edges)
            ft.run_task({})
            self.assertEqual(len(FWAction_patch.call_args[1]["additions"]),12)
        
    def test_babel_PC_frag1(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction") as FWAction_patch:
            pc_frag1 = Molecule.from_file(os.path.join(module_dir, "..", "..", "test_files", "PC_frag1.xyz"))
            ft = FragmentMolecule(molecule=pc_frag1)
            ft.run_task({})
            self.assertEqual(len(FWAction_patch.call_args[1]["additions"]),12)

if __name__ == "__main__":
    unittest.main()
