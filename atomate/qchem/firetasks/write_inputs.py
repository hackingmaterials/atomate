# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

# This module defines firetasks for writing QChem input files

import os

from atomate.utils.utils import load_class
from fireworks import FiretaskBase, explicit_serialize
from pymatgen.io.qchem_io.inputs import QCInput

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
        qchem_input_params (dict): When using a string name for QChem input set, use this as a dict
        to specify kwargs for instantiating the input set parameters. For example, if you want
        to change the DFT_rung, you should provide: {"DFT_rung": ...}.
        This setting is ignored if you provide the full object representation of a QChemDictSet
        rather than a String.
        molecule (Molecule):
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
        # these if statements might need to be reordered at some point
        # if a full QChemDictSet object was provided
        if hasattr(self["qchem_input_set"], "write_file"):
            qcin = self["qchem_input_set"]
        # if a molecule is being passed through fw_spec
        elif fw_spec.get("prev_calc_molecule"):
            mol = fw_spec.get("prev_calc_molecule")
            qcin_cls = load_class("pymatgen.io.qchem_io.sets",
                                  self["qchem_input_set"])
            qcin = qcin_cls(mol, **self.get("qchem_input_params", {}))
        # if a molecule is included as an optional parameter
        elif self.get("molecule"):
            qcin_cls = load_class("pymatgen.io.qchem_io.sets",
                                  self["qchem_input_set"])
            qcin = qcin_cls(
                self.get("molecule"), **self.get("qchem_input_params", {}))
        # if no molecule is present raise an error
        else:
            raise KeyError(
                "No molecule present, add as an optional param or check fw_spec"
            )
        qcin.write_file(input_file)


@explicit_serialize
class WriteCustomInput(FiretaskBase):
    """
        Writes QChem Input files from custom input sets. This firetask gives the maximum flexibility when trying
        to define custom input parameters.

        required_params:
            qchem_input_custom (dict): Define custom input parameters to generate a qchem input file.
            This should be a dictionary of dictionaries (i.e. {{"rem": {"method": "b3lyp", basis": "6-31*G++", ...}
            Each QChem section should be a key with its own dictionary as the value. For more details on how
            the input should be structured look at pymatgen.io.qchem_io.inputs
            ***  ***

        optional_params:
            molecule (Molecule):
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
        # these if statements might need to be reordered at some point
        if "molecule" in self:
            molecule = self["molecule"]
        elif fw_spec.get("prev_calc_molecule"):
            molecule = fw_spec.get("prev_calc_molecule")
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
            molecule=molecule,
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
