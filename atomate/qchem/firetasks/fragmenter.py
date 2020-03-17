# coding: utf-8


# This module defines a task that returns all fragments of a molecule

import copy
from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.analysis.fragmenter import Fragmenter
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
__credits__ = "John Dagdelen, Shyam Dwaraknath, Evan Spotte-Smith"


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
    Thus, we want to simulate charges 0 and 1 of each fragment, given charge=1. Generalizing to
    any positively charged principle with charge P, we simulate P and P-1.

    Realistic fragmentation of a negatively charged molecule (using charge=-1 as an example here)
    will almost always be one of the following:
      molecule(charge=-1) -> fragment1(charge=0) + fragment2(charge=-1)
      molecule(charge=-1) -> fragment1(charge=-1) + fragment2(charge=0)
    Thus, we want to simulate charges -1 and 0 of each fragment, given charge=-1. Generalizing to
    any positively charged principle with charge P, we simulate P and P+1.

    If additional charges are desired by the user, they can be specified with the additional_charges
    input parameter as described below.


    Optional params:
        molecule (Molecule): The molecule to fragment
        edges (list): List of index pairs that define graph edges, aka molecule bonds. If not
                      set, edges will be determined with OpenBabel.
        depth (int): The number of levels of iterative fragmentation to perform, where each
                     level will include fragments obtained by breaking one bond of a fragment
                     one level up. Defaults to 1. However, if set to 0, instead all possible
                     fragments are generated using an alternative, non-iterative scheme.
        open_rings (bool): Whether or not to open any rings encountered during fragmentation.
                           Defaults to True. If true, any bond that fails to yield disconnected
                           graphs when broken is instead removed and the entire structure is
                           optimized with OpenBabel in order to obtain a good initial guess for
                           an opened geometry that can then be put back into QChem to be
                           optimized without the ring just reforming.
        opt_steps (int): Number of optimization steps when opening rings. Defaults to 10000.
        additional_charges (list): List of additional charges besides the defaults described
                                   above. For example, if a principle molecule with a +2 charge
                                   is provided, by default all fragments will be calculated with
                                   +1 and +2 charges as explained above. If the user includes
                                   additional_charges=[0] then all fragments will be calculated
                                   with 0, +1, and +2 charges. Additional charge values of 1 or 2
                                   would not cause any new charges to be calculated as they are
                                   already done. Defaults to [].
        do_triplets (bool): Whether to simulate triplets as well as singlets for molecules with
                            an even number of electrons. Defaults to False.
        linked (bool): If True, use a linked FFopt. If False, use the original. Defaults to True.
        qchem_input_params (dict): Specify kwargs for instantiating the input set parameters.
                                   Basic uses would be to modify the default inputs of the set,
                                   such as dft_rung, basis_set, pcm_dielectric, scf_algorithm,
                                   or max_scf_cycles. See pymatgen/io/qchem/sets.py for default
                                   values of all input parameters. For instance, if a user wanted
                                   to use a more advanced DFT functional, include a pcm with a
                                   dielectric of 30, and use a larger basis, the user would set
                                   qchem_input_params = {"dft_rung": 5, "pcm_dielectric": 30,
                                   "basis_set": "6-311++g**"}. However, more advanced customization
                                   of the input is also possible through the overwrite_inputs key
                                   which allows the user to directly modify the rem, pcm, smd, and
                                   solvent dictionaries that QChemDictSet passes to inputs.py to
                                   print an actual input file. For instance, if a user wanted to
                                   set the sym_ignore flag in the rem section of the input file
                                   to true, then they would set qchem_input_params = {"overwrite_inputs":
                                   "rem": {"sym_ignore": "true"}}. Of course, overwrite_inputs
                                   could be used in conjuction with more typical modifications,
                                   as seen in the test_double_FF_opt workflow test.
        db_file (str): Path to file containing the database credentials. Supports env_chk.
        check_db (bool): Whether or not to check if fragments are present in the database.
                         Defaults to bool(db_file), aka true if a db_file is present and
                         false if db_file is None.
    """

    optional_params = [
        "molecule", "edges", "depth", "open_rings", "opt_steps", "additional_charges", "do_triplets", "linked", "qchem_input_params", "db_file", "check_db"
    ]

    def run_task(self, fw_spec):
        # if a molecule is being passed through fw_spec
        if fw_spec.get("prev_calc_molecule"):
            molecule = fw_spec.get("prev_calc_molecule")
        # if a molecule is included as an optional parameter
        elif self.get("molecule"):
            molecule = self.get("molecule")
        # if no molecule is present raise an error
        else:
            raise KeyError(
                "No molecule present, add as an optional param or check fw_spec"
            )

        self.depth = self.get("depth", 1)
        additional_charges = self.get("additional_charges", [])
        self.do_triplets = self.get("do_triplets", False)
        self.linked = self.get("linked", True)
        self.qchem_input_params = self.get("qchem_input_params", {})

        # Specify charges to consider based on charge of the principle molecule:
        if molecule.charge == 0:
            self.charges = [-1, 0, 1]
        elif molecule.charge > 0:
            self.charges = [molecule.charge-1, molecule.charge]
        else:
            self.charges = [molecule.charge, molecule.charge+1]
        self.principle_charge = molecule.charge

        # Include any additional charges specified by the user:
        for additional_charge in additional_charges:
            if additional_charge not in self.charges:
                print("Adding additional charge " + str(additional_charge))
                self.charges.append(additional_charge)
            else:
                print("Charge " + str(additional_charge) + " already present!")

        # Obtain fragments from Pymatgen's fragmenter:
        fragmenter = Fragmenter(molecule=molecule, edges=self.get("edges", None), depth=self.depth, open_rings=self.get("open_rings", True), opt_steps=self.get("opt_steps", 10000))
        self.unique_fragments = []
        for key in fragmenter.unique_frag_dict:
            for frag in fragmenter.unique_frag_dict[key]:
                self.unique_fragments.append(frag)

        # Convert fragment molecule graphs into molecule objects with charges given in self.charges
        self._build_unique_relevant_molecules()

        # Then find all unique formulae in our unique molecules to facilitate easier database searching
        self.unique_formulae = []
        for molecule in self.unique_molecules:
            if molecule.composition.reduced_formula not in self.unique_formulae:
                self.unique_formulae.append(molecule.composition.reduced_formula)

        # attempt to connect to the database to later check if a fragment has already been calculated
        find_dict = {"formula_pretty": {"$in": self.unique_formulae}}
        if "pcm_dielectric" in self.qchem_input_params:
            find_dict["calcs_reversed.input.solvent.dielectric"] = str(self.qchem_input_params["pcm_dielectric"])
        db_file = env_chk(self.get("db_file"), fw_spec)
        self.check_db = self.get("check_db", bool(db_file))
        self.all_relevant_docs = []
        if db_file and self.check_db:
            mmdb = QChemCalcDb.from_db_file(db_file, admin=True)
            self.all_relevant_docs = list(
                mmdb.collection.find(find_dict, {
                    "formula_pretty": 1,
                    "input.initial_molecule": 1
                }))

        # Return an FWAction which includes a new additional firework for each unique, relevant molecule
        # not already present in our database
        return FWAction(additions=self._build_new_FWs())

    def _build_unique_relevant_molecules(self):
        """
        Convert unique fragments, aka molecule graph objects, into molecule objects with
        a range of charges as defined by self.charges. Builds self.unique_molecules, which
        is then used to generate the new fireworks.
        """
        self.unique_molecules = []
        for unique_fragment in self.unique_fragments:
            for charge in self.charges:
                this_molecule = copy.deepcopy(unique_fragment.molecule)
                this_molecule.set_charge_and_spin(charge=charge)
                self.unique_molecules.append(this_molecule)
        if self.do_triplets:
            for unique_molecule in self.unique_molecules:
                if unique_molecule.spin_multiplicity == 1:
                    this_molecule = copy.deepcopy(unique_molecule)
                    this_molecule.set_charge_and_spin(charge=this_molecule.charge, spin_multiplicity=3)
                    self.unique_molecules.append(this_molecule)

    def _in_database(self, molecule):
        """
        Check if a molecule is already present in the database, which has already been
        queried on relevant formulae and narrowed to self.all_relevant_docs.
        If no docs present, assume fragment is not present
        """
        if len(self.all_relevant_docs) == 0:
            return False

        # otherwise, look through the docs for an entry with an isomorphic molecule with
        # equivalent charge and multiplicity
        else:
            new_mol_graph = MoleculeGraph.with_local_env_strategy(molecule, OpenBabelNN())
            for doc in self.all_relevant_docs:
                if molecule.composition.reduced_formula == doc["formula_pretty"]:
                    old_mol = Molecule.from_dict(doc["input"]["initial_molecule"])
                    old_mol_graph = MoleculeGraph.with_local_env_strategy(old_mol, OpenBabelNN())
                    # If such an equivalent molecule is found, return true
                    if new_mol_graph.isomorphic_to(old_mol_graph) and molecule.charge == old_mol_graph.molecule.charge and molecule.spin_multiplicity == old_mol_graph.molecule.spin_multiplicity:
                        return True
            # Otherwise, return false
            return False

    def _build_new_FWs(self):
        """
        Build the list of new fireworks: a FrequencyFlatteningOptimizeFW for each unique fragment
        molecule, unless the fragment is a single atom, in which case add a SinglePointFW instead.
        If the fragment is already in the database, don't add any new firework.
        """
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
                            max_cores=">>max_cores<<",
                            qchem_input_params=self.qchem_input_params,
                            db_file=">>db_file<<"))
                else:
                    new_FWs.append(
                        FrequencyFlatteningOptimizeFW(
                            molecule=unique_molecule,
                            name="fragment_" + str(ii),
                            qchem_cmd=">>qchem_cmd<<",
                            max_cores=">>max_cores<<",
                            qchem_input_params=self.qchem_input_params,
                            linked=self.linked,
                            db_file=">>db_file<<"))
        return new_FWs
