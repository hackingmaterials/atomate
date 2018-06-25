# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

# This module defines a task that returns all fragments of a molecule


from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from atomate.utils.utils import env_chk
from itertools import combinations
import networkx as nx
from fireworks import FiretaskBase, FWAction, explicit_serialize

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


@explicit_serialize
class FragmentMolecule(FiretaskBase):
    """
    Find all unique fragments of a molecule

    Optional params:
        molecule (Molecule): 
        edges (list): List of index pairs that define graph edges, aka molecule bonds
        max_cores (int): Maximum number of cores to parallelize over. Defaults to 32.
        qchem_input_params (dict): Specify kwargs for instantiating the input set parameters.
                                   For example, if you want to change the DFT_rung, you should
                                   provide: {"DFT_rung": ...}. Defaults to None.
        db_file (str): path to file containing the database credentials. Supports env_chk.

    """

    optional_params = ["molecule", "edges", "max_cores", "qchem_input_params", "db_file"]

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
        
        # if edges are not passed by the user and babel is not available, raise an error
        if not self.get("edges") and not have_babel:
            raise KeyError(
                "OpenBabel not accessible and no bonds provided! Exiting..."
            )

        # build the MoleculeGraph
        mol_graph = build_MoleculeGraph(mol, self.get("edges", None))

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
        from atomate.qchem.fireworks.core import FrequencyFlatteningOptimizeFW
        from atomate.qchem.fireworks.core import SinglePointFW
        new_FWs = []
        for ii,unique_molecule in enumerate(unique_molecules):
            if not_in_database(unique_molecule):
                if len(unique_molecule) == 1:
                    new_FWs.append(SinglePointFW(molecule=unique_molecule,
                                                 name="fragment_"+str(ii),
                                                 qchem_cmd=">>qchem_cmd<<",
                                                 max_cores=self.get("max_cores", 32),
                                                 qchem_input_params=self.get("qchem_input_params", {}),
                                                 db_file=">>db_file<<"))
                else:
                    new_FWs.append(FrequencyFlatteningOptimizeFW(molecule=unique_molecule,
                                                                 name="fragment_"+str(ii),
                                                                 qchem_cmd=">>qchem_cmd<<",
                                                                 max_cores=self.get("max_cores", 32),
                                                                 qchem_input_params=self.get("qchem_input_params", {}),
                                                                 db_file=">>db_file<<"))

        return FWAction(additions=new_FWs)
        

def build_MoleculeGraph(molecule, edges):
    if edges == None:
        babel_mol = BabelMolAdaptor(molecule).openbabel_mol
        edges = []
        for obbond in ob.OBMolBondIter(babel_mol):
            edges += [[obbond.GetBeginAtomIdx()-1,obbond.GetEndAtomIdx()-1]]
    mol_graph = MoleculeGraph.with_empty_graph(molecule)
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
    return mol_graph

def _node_match(node, othernode):
    return node["specie"] == othernode["specie"]

def is_isomorphic(graph1, graph2):
    return nx.is_isomorphic(graph1, graph2, node_match=_node_match)

def not_in_database(molecule):
    # get the database connection info
    db_file = env_chk(self.get("db_file"), fw_spec)

    # if we cannot connect to the database, assume fragment is not present
    if not db_file:
        return True

    # otherwise, connect to the database
    else:
        mmdb = QChemCalcDb.from_db_file(db_file, admin=True)

        new_mol_graph = build_MoleculeGraph(molecule, None)
        for doc in mmdb.collection.find({"formula_pretty": molecule.composition.reduced_formula}):
            old_mol_graph = build_MoleculeGraph(doc["output"]["initial_molecule"], None)
            if nx.is_isomorphic(new_mol_graph.graph, old_mol_graph.graph):
                return False
        return True

    