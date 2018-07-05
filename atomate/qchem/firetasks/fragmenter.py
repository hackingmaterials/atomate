# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

# This module defines a task that returns all fragments of a molecule

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from atomate.utils.utils import env_chk
from atomate.qchem.database import QChemCalcDb
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

    optional_params = [
        "molecule", "edges", "max_cores", "qchem_input_params", "db_file"
    ]

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
                "OpenBabel not accessible and no bonds provided! Exiting...")

        # build the MoleculeGraph
        mol_graph = build_MoleculeGraph(mol, self.get("edges", None))

        # find all unique fragments
        unique_fragments = build_unique_fragments(mol_graph)

        # build three molecule objects for each unique fragment:
        # original charge, original charge +1, original charge -1
        unique_molecules = build_unique_molecules(unique_fragments, mol.charge)

        unique_formulae = []
        for molecule in unique_molecules:
            if molecule.composition.reduced_formula not in unique_formulae:
                unique_formulae.append(molecule.composition.reduced_formula)

        # attempt to connect to the database to later check if a fragment has already been calculated
        db_file = env_chk(self.get("db_file"), fw_spec)
        all_relevant_docs = []
        # if db_file:
        #     mmdb = QChemCalcDb.from_db_file(db_file, admin=True)
        #     all_relevant_docs = list(
        #         mmdb.collection.find({
        #             "formula_pretty": {
        #                 "$in": unique_formulae
        #             }
        #         }, {
        #             "formula_pretty": 1,
        #             "output.initial_molecule": 1
        #         }))

        # build the list of new fireworks
        new_FWs = build_new_FWs(unique_molecules, all_relevant_docs,
                                self.get("max_cores", 32),
                                self.get("qchem_input_params", {}))

        return FWAction(additions=new_FWs)


def edges_from_babel(molecule):
    babel_mol = BabelMolAdaptor(molecule).openbabel_mol
    edges = []
    for obbond in ob.OBMolBondIter(babel_mol):
        edges += [[obbond.GetBeginAtomIdx() - 1, obbond.GetEndAtomIdx() - 1]]
    return edges


def build_MoleculeGraph(molecule, edges=None):
    if edges == None:
        edges = edges_from_babel(molecule)
    mol_graph = MoleculeGraph.with_empty_graph(molecule)
    for edge in edges:
        mol_graph.add_edge(edge[0], edge[1])
    mol_graph.graph = mol_graph.graph.to_undirected()
    species = {}
    coords = {}
    for node in mol_graph.graph:
        species[node] = mol_graph.molecule[node].specie.symbol
        coords[node] = mol_graph.molecule[node].coords
    nx.set_node_attributes(mol_graph.graph, species, "specie")
    nx.set_node_attributes(mol_graph.graph, coords, "coords")
    return mol_graph


def build_unique_fragments(mol_graph):
    # find all possible fragments, aka connected induced subgraphs
    all_fragments = []
    for ii in range(1, len(mol_graph.molecule)):
        for combination in combinations(mol_graph.graph.nodes, ii):
            subgraph = nx.subgraph(mol_graph.graph, combination)
            if nx.is_connected(subgraph):
                all_fragments.append(subgraph)

    # narrow to all unique fragments using graph isomorphism
    unique_fragments = []
    for fragment in all_fragments:
        if not [is_isomorphic(fragment, f)
                for f in unique_fragments].count(True) >= 1:
            unique_fragments.append(fragment)
    return unique_fragments


def build_unique_molecules(unique_fragments, orig_charge):
    unique_molecules = []
    for fragment in unique_fragments:
        species = [fragment.node[ii]["specie"] for ii in fragment.nodes]
        coords = [fragment.node[ii]["coords"] for ii in fragment.nodes]
        unique_molecule0 = Molecule(
            species=species, coords=coords, charge=orig_charge)
        unique_molecule1 = Molecule(
            species=species, coords=coords, charge=orig_charge + 1)
        unique_molecule2 = Molecule(
            species=species, coords=coords, charge=orig_charge - 1)
        unique_molecules.append(unique_molecule0)
        unique_molecules.append(unique_molecule1)
        unique_molecules.append(unique_molecule2)
    return unique_molecules


def _node_match(node, othernode):
    return node["specie"] == othernode["specie"]


def is_isomorphic(graph1, graph2):
    return nx.is_isomorphic(graph1, graph2, node_match=_node_match)


def not_in_database(molecule, docs):
    # if no docs present, assume fragment is not present
    if len(docs) == 0:
        return True

    # otherwise, look through the docs for an equivalent entry
    else:
        new_mol_graph = build_MoleculeGraph(molecule, None)
        for doc in docs:
            if molecule.composition.reduced_formula == doc["formula_pretty"]:
                try:
                    old_mol = Molecule.from_dict(doc["output"]["initial_molecule"])
                except TypeError:
                    old_mol = doc["output"]["initial_molecule"]
                old_mol_graph = build_MoleculeGraph(old_mol,None)
                if nx.is_isomorphic(
                        new_mol_graph.graph, old_mol_graph.graph
                ) and molecule.charge == old_mol_graph.molecule.charge and molecule.spin_multiplicity == old_mol_graph.molecule.spin_multiplicity:
                    return False
        return True


def build_new_FWs(unique_molecules, all_relevant_docs, max_cores,
                  qchem_input_params):
    # build the list of new fireworks: a FrequencyFlatteningOptimizeFW for each unique fragment
    # unless the fragment is a single atom, in which case add a SinglePointFW instead.
    from atomate.qchem.fireworks.core import FrequencyFlatteningOptimizeFW
    from atomate.qchem.fireworks.core import SinglePointFW
    new_FWs = []
    for ii, unique_molecule in enumerate(unique_molecules):
        if not_in_database(unique_molecule, all_relevant_docs):
            if len(unique_molecule) == 1:
                new_FWs.append(
                    SinglePointFW(
                        molecule=unique_molecule,
                        name="fragment_" + str(ii),
                        qchem_cmd=">>qchem_cmd<<",
                        max_cores=max_cores,
                        qchem_input_params=qchem_input_params,
                        db_file=">>db_file<<"))
            else:
                new_FWs.append(
                    FrequencyFlatteningOptimizeFW(
                        molecule=unique_molecule,
                        name="fragment_" + str(ii),
                        qchem_cmd=">>qchem_cmd<<",
                        max_cores=max_cores,
                        qchem_input_params=qchem_input_params,
                        db_file=">>db_file<<"))
    return new_FWs
