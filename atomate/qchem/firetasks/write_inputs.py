# coding: utf-8


# This module defines firetasks for writing QChem input files

import os

from atomate.utils.utils import load_class
from fireworks import FiretaskBase, explicit_serialize
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN

__author__ = "Brandon Wood"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Brandon Wood"
__email__ = "b.wood@berkeley.edu"
__status__ = "Alpha"
__date__ = "5/20/18"
__credits__ = "Sam Blau, Shyam Dwaraknath"


@explicit_serialize
class WriteInputFromIOSet(FiretaskBase):
    """
    Writes QChem Input files from input sets. A dictionary is passed to WriteInputFromIOSet where
    parameters are given as keys in the dictionary.

    required_params:
        qc_input_set (QChemDictSet or str): Either a QChemDictSet object or a string
        name for the QChem input set (e.g., "OptSet"). *** Note that if the molecule is to be inherited through
        fw_spec qc_input_set must be a string name for the QChem input set. ***

    optional_params:
        molecule (Molecule): Molecule that will be subjected to an electronic structure calculation
        qchem_input_params (dict): When using a string name for QChem input set, use this as a dict
                                   to specify kwargs for instantiating the input set parameters. This
                                   setting is ignored if you provide the full object representation of
                                   a QChemDictSet. Basic uses would be to modify the default inputs of
                                   the set, such as dft_rung, basis_set, pcm_dielectric, scf_algorithm,
                                   or max_scf_cycles. See pymatgen/io/qchem/sets.py for default values
                                   of all input parameters. For instance, if a user wanted to use a more
                                   advanced DFT functional, include a pcm with a dielectric of 30, and
                                   use a larger basis, the user would set qchem_input_params = {"dft_rung":
                                   5, "pcm_dielectric": 30, "basis_set": "6-311++g**"}. However, more
                                   advanced customization of the input is also possible through the
                                   overwrite_inputs key which allows the user to directly modify the rem,
                                   pcm, smd, and solvent dictionaries that QChemDictSet passes to inputs.py
                                   to print an actual input file. For instance, if a user wanted to set the
                                   sym_ignore flag in the rem section of the input file to true, then they
                                   would set qchem_input_params = {"overwrite_inputs": "rem": {"sym_ignore":
                                   "true"}}. Of course, overwrite_inputs could be used in conjuction with
                                   more typical modifications, as seen in the test_double_FF_opt workflow
                                   test.
        input_file (str): Name of the QChem input file. Defaults to mol.qin
        write_to_dir (str): Path of the directory where the QChem input file will be written,
        the default is to write to the current working directory
    """

    required_params = ["qchem_input_set"]
    optional_params = [
        "molecule", "qchem_input_params", "input_file", "write_to_dir"
    ]

    def run_task(self, fw_spec):
        input_file = os.path.join(self.get("write_to_dir", ""),self.get("input_file", "mol.qin"))

        # if a full QChemDictSet object was provided
        if hasattr(self["qchem_input_set"], "write_file"):
            qcin = self["qchem_input_set"]
        # if a molecule is being passed through fw_spec
        elif fw_spec.get("prev_calc_molecule"):
            prev_calc_mol = fw_spec.get("prev_calc_molecule")
            # if a molecule is also passed as an optional parameter
            if self.get("molecule"):
                mol = self.get("molecule")
                # check if mol and prev_calc_mol are isomorphic
                mol_graph = MoleculeGraph.with_local_env_strategy(mol,
                                                                  OpenBabelNN(),
                                                                  reorder=False,
                                                                  extend_structure=False)
                prev_mol_graph = MoleculeGraph.with_local_env_strategy(prev_calc_molecule,
                                                                       OpenBabelNN(),
                                                                       reorder=False,
                                                                       extend_structure=False)
                # If they are isomorphic, aka a previous FW has not changed bonding,
                # then we will use prev_calc_mol. If bonding has changed, we will use mol.
                if mol_graph.isomorphic_to(prev_mol_graph):
                    mol = prev_calc_mol
                elif self["qchem_input_set"] != "OptSet":
                    print("WARNING: Molecule from spec is not isomorphic to passed molecule!")
                    mol = prev_calc_mol
                else:
                    print("Not using prev_calc_mol as it is not isomorphic to passed molecule!")
            else:
              mol = prev_calc_mol

            qcin_cls = load_class("pymatgen.io.qchem.sets",
                                  self["qchem_input_set"])
            qcin = qcin_cls(mol, **self.get("qchem_input_params", {}))
        # if a molecule is only included as an optional parameter
        elif self.get("molecule"):
            qcin_cls = load_class("pymatgen.io.qchem.sets",
                                  self["qchem_input_set"])
            qcin = qcin_cls(
                self.get("molecule"), **self.get("qchem_input_params", {}))
        # if no molecule is present raise an error
        else:
            raise KeyError(
                "No molecule present, add as an optional param or check fw_spec"
            )
        qcin.write(input_file)


@explicit_serialize
class WriteCustomInput(FiretaskBase):
    """
        Writes QChem Input files from custom input sets. This firetask gives the maximum flexibility when trying
        to define custom input parameters.

        required_params:
            qchem_input_custom (dict): Define custom input parameters to generate a qchem input file.
            This should be a dictionary of dictionaries (i.e. {{"rem": {"method": "b3lyp", basis": "6-31*G++", ...}
            Each QChem section should be a key with its own dictionary as the value. For more details on how
            the input should be structured look at pymatgen.io.qchem.inputs
            ***  ***

        optional_params:
            input_file (str): Name of the QChem input file. Defaults to mol.qin
            write_to_dir (str): Path of the directory where the QChem input file will be written,
            the default is to write to the current working directory
        """

    required_params = ["rem"]
    # optional_params will need to be modified if more QChem sections are added QCInput
    optional_params = [
        "molecule", "opt", "pcm", "solvent", "input_file", "write_to_dir"
    ]

    def run_task(self, fw_spec):
        input_file = os.path.join(self.get("write_to_dir", ""),self.get("input_file", "mol.qin"))
        # if a molecule is being passed through fw_spec
        if fw_spec.get("prev_calc_molecule"):
            prev_calc_mol = fw_spec.get("prev_calc_molecule")
            # if a molecule is also passed as an optional parameter
            if self.get("molecule"):
                mol = self.get("molecule")
                # check if mol and prev_calc_mol are isomorphic
                mol_graph = MoleculeGraph.with_local_env_strategy(mol,
                                                                  OpenBabelNN(),
                                                                  reorder=False,
                                                                  extend_structure=False)
                prev_mol_graph = MoleculeGraph.with_local_env_strategy(prev_calc_molecule,
                                                                       OpenBabelNN(),
                                                                       reorder=False,
                                                                       extend_structure=False)
                if mol_graph.isomorphic_to(prev_mol_graph):
                    mol = prev_calc_mol
                else:
                    print("WARNING: Molecule from spec is not isomorphic to passed molecule!")
            else:
              mol = prev_calc_mol
        elif self.get("molecule"):
            mol = self.get("molecule")
        else:
            raise KeyError(
                "No molecule present, add as an optional param or check fw_spec"
            )
        # in the current structure there needs to be a statement for every optional QChem section
        # the code below defaults the section to None if the variable is not passed
        opt = self.get("opt", None)
        pcm = self.get("pcm", None)
        solvent = self.get("solvent", None)

        qcin = QCInput(
            molecule=mol,
            rem=self["rem"],
            opt=opt,
            pcm=pcm,
            solvent=solvent)
        qcin.write_file(input_file)


@explicit_serialize
class WriteInput(FiretaskBase):
    """
    Writes QChem input file from QCInput object.

    required_params:
        qc_input (QCInput): QCInput object

    optional_params:
        input_file (str): Name of the QChem input file. Defaults to mol.qin
        write_to_dir (str): Path of the directory where the QChem input file will be written,
        the default is to write to the current working directory

    """
    required_params = ["qc_input"]
    optional_params = ["input_file", "write_to_dir"]

    def run_task(self, fw_spec):
        # if a QCInput object is provided
        input_file = os.path.join(self.get("write_to_dir", ""),self.get("input_file", "mol.qin"))

        qcin = self["qc_input"]
        qcin.write_file(input_file)
