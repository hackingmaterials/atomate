# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

# This module defines a task that returns all fragments of a molecule

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import build_MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from atomate.utils.utils import env_chk
from atomate.qchem.database import QChemCalcDb
from itertools import combinations
import networkx as nx
from fireworks import FiretaskBase, FWAction, explicit_serialize

have_babel = True
try:
    from pymatgen.io.babel import BabelMolAdaptor
    import openbabel as ob
except ImportError:
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
        "molecule", "edges", "max_cores", "qchem_input_params", "db_file", "check_db"
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
        edges = self.get("edges", None)
        if edges is None:
            mol_graph = build_MoleculeGraph(mol, strategy=OpenBabelNN,
                                            reorder=False, extend_structure=False)
        else:
            mol_graph = build_MoleculeGraph(mol, edges=edges)

        # find all unique fragments
        unique_fragments = mol_graph.build_unique_fragments()

        # build three molecule objects for each unique fragment:
        # original charge, original charge +1, original charge -1
        unique_molecules = build_unique_molecules(unique_fragments, mol.charge)

        unique_formulae = []
        for molecule in unique_molecules:
            if molecule.composition.reduced_formula not in unique_formulae:
                unique_formulae.append(molecule.composition.reduced_formula)

        # attempt to connect to the database to later check if a fragment has already been calculated
        db_file = env_chk(self.get("db_file"), fw_spec)
        check_db = self.get("check_db", True)
        all_relevant_docs = []
        if db_file and check_db:
            mmdb = QChemCalcDb.from_db_file(db_file, admin=True)
            all_relevant_docs = list(
                mmdb.collection.find({
                    "formula_pretty": {
                        "$in": unique_formulae
                    }
                }, {
                    "formula_pretty": 1,
                    "output.initial_molecule": 1
                }))

        # build the list of new fireworks
        new_FWs = build_new_FWs(unique_molecules, all_relevant_docs,
                                self.get("max_cores", 32),
                                self.get("qchem_input_params", {}))

        return FWAction(additions=new_FWs)


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


def not_in_database(molecule, docs):
    # if no docs present, assume fragment is not present
    if len(docs) == 0:
        return True

    # otherwise, look through the docs for an equivalent entry
    else:
        new_mol_graph = build_MoleculeGraph(molecule, strategy=OpenBabelNN,
                                            reorder=False, extend_structure=False)
        for doc in docs:
            if molecule.composition.reduced_formula == doc["formula_pretty"]:
                try:
                    old_mol = Molecule.from_dict(doc["output"]["initial_molecule"])
                except TypeError:
                    old_mol = doc["output"]["initial_molecule"]
                old_mol_graph = build_MoleculeGraph(old_mol, strategy=OpenBabelNN,
                                                    reorder=False, extend_structure=False)
                if new_mol_graph.isomorphic_to(old_mol_graph) and molecule.charge == old_mol_graph.molecule.charge and molecule.spin_multiplicity == old_mol_graph.molecule.spin_multiplicity:
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
