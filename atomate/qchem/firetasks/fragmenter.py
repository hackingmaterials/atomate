# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

# This module defines a task that returns all fragments of a molecule

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import build_MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from atomate.utils.utils import env_chk
from atomate.qchem.database import QChemCalcDb
from fireworks import FiretaskBase, FWAction, explicit_serialize

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
    Simulate all unique fragments of a molecule relevant for realistic fragmentation
    of the molecule itself or of its subfragments.

    Realistic fragmentation of a neutral molecule will almost always be one of the following:
      molecule(charge=0) -> fragment1(charge=0) + fragment2(charge=0)
      molecule(charge=0) -> fragment1(charge=1) + fragment2(charge=-1)
      molecule(charge=0) -> fragment1(charge=-1) + fragment2(charge=1)
    Thus, we want to simulate charges -1, 0, and 1 of each fragment.

    Realistic fragmentation of a positively charged molecule (using charge=1 as an example here)
    will almost always be one of the following:
      molecule(charge=1) -> fragment1(charge=0) + fragment2(charge=1)
      molecule(charge=1) -> fragment1(charge=1) + fragment2(charge=0)
      molecule(charge=1) -> fragment1(charge=2) + fragment2(charge=-1)
      molecule(charge=1) -> fragment1(charge=-1) + fragment2(charge=2)
    where the last two are significantly less likely, but possible.
    Thus, we want to simulate charges -1, 0, 1, and 2 of each fragment, given charge=1.
    We generalize to any positive charge in _build_unique_relevant_molecules.
    
    Realistic fragmentation of a negatively charged molecule (using charge=-1 as an example here)
    will almost always be one of the following:
      molecule(charge=-1) -> fragment1(charge=0) + fragment2(charge=-1)
      molecule(charge=-1) -> fragment1(charge=-1) + fragment2(charge=0)
      molecule(charge=-1) -> fragment1(charge=-2) + fragment2(charge=1)
      molecule(charge=-1) -> fragment1(charge=1) + fragment2(charge=-2)
    where the last two are significantly less likely, but possible.
    Thus, we want to simulate charges -2, -1, 0, and 1 of each fragment, given charge=-1.
    We generalize to any negative charge in _build_unique_relevant_molecules.
    

    Optional params:
        molecule (Molecule): The molecule to fragment
        edges (list): List of index pairs that define graph edges, aka molecule bonds
        max_cores (int): Maximum number of cores to parallelize over. Defaults to 32.
        qchem_input_params (dict): Specify kwargs for instantiating the input set parameters.
                                   For example, if you want to change the DFT_rung, you should
                                   provide: {"DFT_rung": ...}. Defaults to None.
        db_file (str): Path to file containing the database credentials. Supports env_chk.
        check_db (bool): Whether or not to check if fragments are present in the database.
                         Defaults to bool(db_file), aka true if a db_file is present.
    """

    optional_params = [
        "molecule", "edges", "max_cores", "qchem_input_params", "db_file", "check_db"
    ]

    def run_task(self, fw_spec):
        # if a molecule is being passed through fw_spec
        if fw_spec.get("prev_calc_molecule"):
            self.mol = fw_spec.get("prev_calc_molecule")
        # if a molecule is included as an optional parameter
        elif self.get("molecule"):
            self.mol = self.get("molecule")
        # if no molecule is present raise an error
        else:
            raise KeyError(
                "No molecule present, add as an optional param or check fw_spec"
            )

        # build the MoleculeGraph
        edges = self.get("edges", None)
        if edges is None:
            mol_graph = build_MoleculeGraph(self.mol, strategy=OpenBabelNN,
                                            reorder=False, extend_structure=False)
        else:
            mol_graph = build_MoleculeGraph(self.mol, edges=edges)

        # Find all unique fragments
        self.unique_fragments = mol_graph.build_unique_fragments()

        # Build all unique molecules relevant to realistic fragmentation
        self._build_unique_relevant_molecules()

        self.unique_formulae = []
        for molecule in self.unique_molecules:
            if molecule.composition.reduced_formula not in self.unique_formulae:
                self.unique_formulae.append(molecule.composition.reduced_formula)

        # attempt to connect to the database to later check if a fragment has already been calculated
        db_file = env_chk(self.get("db_file"), fw_spec)
        check_db = self.get("check_db", bool(db_file))
        self.all_relevant_docs = []
        if db_file and check_db:
            mmdb = QChemCalcDb.from_db_file(db_file, admin=True)
            self.all_relevant_docs = list(
                mmdb.collection.find({
                    "formula_pretty": {
                        "$in": self.unique_formulae
                    }
                }, {
                    "formula_pretty": 1,
                    "output.initial_molecule": 1
                }))

        # Return an FWAction which includes a new additional firework for each unique, relevant molecule
        # not already present in our database
        return FWAction(additions=self._build_new_FWs())

    def _build_unique_relevant_molecules(self):
        # Convert unique fragment objects to unique molecule objects, where each fragment
        # yields three or four molecules given the range of relevant charges for dissociation
        # as explained at the top of this firetask.
        self.unique_molecules = []
        if self.mol.charge == 0:
            charges = [-1, 0, 1]
        elif self.mol.charge > 0:
            charges = [self.mol.charge-2, self.mol.charge-1, self.mol.charge, self.mol.charge+1]
        else:
            charges = [self.mol.charge-1, self.mol.charge, self.mol.charge+1, self.mol.charge+2]
        for fragment in self.unique_fragments:
            species = [fragment.node[ii]["specie"] for ii in fragment.nodes]
            coords = [fragment.node[ii]["coords"] for ii in fragment.nodes]
            for charge in charges:
                self.unique_molecules.append(Molecule(species=species, coords=coords, charge=charge))

    def _in_database(self, molecule):
        # Check if a molecule is already present in the database, which has already been
        # queried on relevant formulae and narrowed to self.all_relevant_docs.
        # If no docs present, assume fragment is not present
        if len(self.all_relevant_docs) == 0:
            return False

        # otherwise, look through the docs for an entry with an isomorphic molecule with
        # equivalent charge and multiplicity
        else:
            new_mol_graph = build_MoleculeGraph(molecule, strategy=OpenBabelNN,
                                                reorder=False, extend_structure=False)
            for doc in self.all_relevant_docs:
                if molecule.composition.reduced_formula == doc["formula_pretty"]:
                    try:
                        old_mol = Molecule.from_dict(doc["output"]["initial_molecule"])
                    except TypeError:
                        old_mol = doc["output"]["initial_molecule"]
                    old_mol_graph = build_MoleculeGraph(old_mol, strategy=OpenBabelNN,
                                                        reorder=False, extend_structure=False)
                    # If such an equivalent molecule is found, return true
                    if new_mol_graph.isomorphic_to(old_mol_graph) and molecule.charge == old_mol_graph.molecule.charge and molecule.spin_multiplicity == old_mol_graph.molecule.spin_multiplicity:
                        return True
            # Otherwise, return false
            return False

    def _build_new_FWs(self):
        # Build the list of new fireworks: a FrequencyFlatteningOptimizeFW for each unique fragment
        # molecule, unless the fragment is a single atom, in which case add a SinglePointFW instead.
        # If the fragment is already in the database, don't add any new firework.
        from atomate.qchem.fireworks.core import FrequencyFlatteningOptimizeFW
        from atomate.qchem.fireworks.core import SinglePointFW
        new_FWs = []
        for ii, unique_molecule in enumerate(self.unique_molecules):
            if not self._in_database(unique_molecule):
                if len(unique_molecule) == 1:
                    new_FWs.append(
                        SinglePointFW(
                            molecule=unique_molecule,
                            name="fragment_" + str(ii),
                            qchem_cmd=">>qchem_cmd<<",
                            max_cores=self.get("max_cores",32),
                            qchem_input_params=self.get("qchem_input_params", {}),
                            db_file=">>db_file<<"))
                else:
                    new_FWs.append(
                        FrequencyFlatteningOptimizeFW(
                            molecule=unique_molecule,
                            name="fragment_" + str(ii),
                            qchem_cmd=">>qchem_cmd<<",
                            max_cores=self.get("max_cores",32),
                            qchem_input_params=self.get("qchem_input_params", {}),
                            db_file=">>db_file<<"))
        return new_FWs
