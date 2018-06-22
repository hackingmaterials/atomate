# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

# This module defines a task that returns all fragments of a molecule


from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from itertools import combinations
import networkx as nx
from fireworks import FiretaskBase, FWAction, explicit_serialize
from atomate.qchem.fireworks.core import FrequencyFlatteningOptimizeFW

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
__credits__ = "John Dagdelen, Shyam Dwaraknath"


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
        
        # if edges are passed by the user
        if self.get("edges"):
            edges = self.get("edges")
        # if not, and if we don't have babel, raise an error
        elif not have_babel:
            raise KeyError(
                "OpenBabel not accessible and no bonds provided! Exiting..."
            )
        # if using babel
        else:
            babel_mol = BabelMolAdaptor(mol).openbabel_mol
            edges = []
            for obbond in ob.OBMolBondIter(babel_mol):
                edges += [[obbond.GetBeginAtomIdx()-1,obbond.GetEndAtomIdx()-1]]

        # build the MoleculeGraph
        mol_graph = MoleculeGraph.with_empty_graph(mol)
        for edge in edges:
            mol_graph.add_edge(edge[0],edge[1])
        mol_graph.graph = mol_graph.graph.to_undirected()
        species = {}
        coords = {}
        for node in mol_graph.graph:
            species[node] = mol_graph.molecule[node].specie.symbol
            coords[node] = mol_graph.molecule[node].coords
        nx.set_node_attributes(mol_graph.graph, species, "specie")
        nx.set_node_attributes(mol_graph.graph, coords, "coords")

        # find all possible fragments, aka connected induced subgraphs
        all_fragments = []
        for ii in range(1,len(mol)):
            for combination in combinations(mol_graph.graph.nodes,ii):
                subgraph = nx.subgraph(mol_graph.graph, combination)
                if nx.is_connected(subgraph):
                    all_fragments.append(subgraph)

        # narrow to all unique fragments using graph isomorphism
        unique_fragments = []
        for fragment in all_fragments:
            if not [is_isomorphic(fragment, f) for f in unique_fragments].count(True) >= 1:
                unique_fragments.append(fragment)
                
        # build three molecule objects for each unique fragment: 
        # original charge, original charge +1, original charge -1
        unique_molecules = []
        for fragment in unique_fragments:
            species = [fragment.node[ii]["specie"] for ii in fragment.nodes]
            coords = [fragment.node[ii]["coords"] for ii in fragment.nodes]
            unique_molecule0 = Molecule(species=species, 
                                       coords=coords, 
                                       charge=mol.charge)
            unique_molecule1 = Molecule(species=species, 
                                       coords=coords, 
                                       charge=mol.charge+1)
            unique_molecule2 = Molecule(species=species, 
                                       coords=coords, 
                                       charge=mol.charge-1)
            unique_molecules.append(unique_molecule0)
            unique_molecules.append(unique_molecule1)
            unique_molecules.append(unique_molecule2)

        # build the list of new fireworks: a FrequencyFlatteningOptimizeFW for each unique fragment
        new_FWs = []
        for unique_molecule in enumerate(unique_molecules):
            new_FWs.append(FrequencyFlatteningOptimizeFW(molecule=unique_molecule))
            
        return FWAction(additions=new_FWs)


def _node_match(node, othernode):
    return node["specie"] == othernode["specie"]

def is_isomorphic(graph1, graph2):
    return nx.is_isomorphic(graph1, graph2, node_match=_node_match)
    