# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

# This module defines a task that returns all fragments of a molecule


from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from fireworks import FiretaskBase, explicit_serialize

have_babel = True
try:
    from pymatgen.io.babel import BabelMolAdaptor
    import openbabel as ob
except:
    print("Cannot find OpenBabel! Thus, bonds must be provided by the user.")
    have_babel = False


__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "6/13/18"


DEBUG = False


@explicit_serialize
class FragmentMolecule(FiretaskBase):
    """
    Find all unique fragments of a molecule

    Optional params:
        molecule (Molecule): 
        edges (list): 

    """

    optional_params = ["molecule", "edges"]

    def run_task(self, fw_spec):
        # if a molecule is being passed through fw_spec
        if fw_spec.get("prev_calc_molecule"):
            mol = fw_spec.get("prev_calc_molecule")
        # if a molecule is included as an optional parameter
        elif self.get("molecule"):
            mol = self.get("molecule")
        # if no molecule is present raise an error
        else:
            raise KeyError(
                "No molecule present, add as an optional param or check fw_spec"
            )
        
        if self.get("edges"):
            edges = self.get("edges")
        elif not have_babel:
            raise KeyError(
                "OpenBabel not accessible and no bonds provided! Exiting..."
            )
        else:
            babel_mol = BabelMolAdaptor(mol).openbabel_mol
            edges = []
            for obbond in ob.OBMolBondIter(babel_mol):
                edges += [[obbond.GetBeginAtomIdx()-1,obbond.GetEndAtomIdx()-1]]

        mol_graph = MoleculeGraph.with_empty_graph(mol)
        for edge in edges:
            mol_graph.add_edge(edge[0],edge[1])

        done = False
        self.all_unique_frags = []

        while not done:
            if len(self.all_unique_frags) == 0:
                frags = fragment(mol_graph)
                new_frags_added = self.add_unique_frags(frags)
            else:
                total_new_frags = []
                for ii,unique_frag in enumerate(self.all_unique_frags):
                    subfrags = fragment(unique_frag)
                    new_frags_added = self.add_unique_frags(subfrags)
                    total_new_frags += new_frags_added
                if len(total_new_frags) == 0:
                    done = True


        # TODO insert into database
        # TODO how am I naming things for easy referencing?

    def add_unique_frags(self, frags):
        new_frags_added = []
        for new_frag in frags:
            is_new = True
            for old_frag in self.all_unique_frags:
                if new_frag.molecule.species == old_frag.molecule.species:
                    if len(new_frag.molecule.species) == 1:
                        is_new = is_new and False
                    else:
                        is_new = is_new and not new_frag.equivalent_to(old_frag)
            if is_new:
                self.all_unique_frags += [new_frag]
                new_frags_added += [new_frag]
        if DEBUG:
            self.check_unique()
        return new_frags_added

    def check_unique(self):
        for ii in range(len(self.all_unique_frags)):
            for jj in range(ii+1, len(self.all_unique_frags)):
                if self.all_unique_frags[ii].molecule.species == self.all_unique_frags[jj].molecule.species:
                    if self.all_unique_frags[ii].equivalent_to(self.all_unique_frags[jj]):
                        print(self.all_unique_frags[ii])
                        print(self.all_unique_frags[jj])
                        raise ValueError("Equivalent fragments found! Exiting...")

def fragment(mol_graph):
    all_frags = []
    ring_bonds = []
    bad_pairs = []
    for bond in mol_graph.graph.edges:
        bonds = [(bond[0],bond[1])]
        try:
            frags = mol_graph.split_molecule_subgraphs(bonds,allow_reverse=True)
            all_frags += frags
        except RuntimeError:
            ring_bonds += bonds

    bond_pairs = []
    for ii,bond in enumerate(ring_bonds):
        for jj in range(ii+1,len(ring_bonds)):
            bond_pair = [bond, ring_bonds[jj]]
            bond_pairs += [bond_pair]

    for bond_pair in bond_pairs:
        try:
            frags = mol_graph.split_molecule_subgraphs(bond_pair,allow_reverse=True)
            all_frags += frags
        except RuntimeError:
            bad_pairs += bond_pair
    return all_frags


