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

from atomate.qchem.firetasks.fragmenter import FragmentMolecule, build_unique_fragments, build_unique_molecules, build_new_FWs, build_MoleculeGraph, is_isomorphic
from atomate.utils.testing import AtomateTest
from pymatgen.core import Molecule
from pymatgen.analysis.graphs import MoleculeGraph

__author__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class TestTop11(AtomateTest):

    def setUp(self, lpad=False):
        super(TestTop11, self).setUp(lpad=False)

    def tearDown(self):
        pass

    def test_top_11(self):
        mol_names = ["BF4-.xyz","DEC.xyz","DMC.xyz","EC.xyz","EMC.xyz","FEC.xyz","FSI-.xyz","PC.xyz","PF6-.xyz","TFSI-.xyz","VC.xyz"]
        all_fragments = []
        for name in mol_names:
            mol = Molecule.from_file(os.path.join(module_dir, "..", "..", "test_files", name))
            mol_graph = build_MoleculeGraph(mol)
            unique_fragments = build_unique_fragments(mol_graph)
            print(len(unique_fragments))
            all_fragments.append(unique_fragments)

        print()
        print(len(all_fragments))
        for fragment in all_fragments:
            if not [is_isomorphic(fragment, f) for f in unique_fragments].count(True) >= 1:
                unique_fragments.append(fragment)
        
        print()
        print(len(unique_fragments))


if __name__ == "__main__":
    unittest.main()
