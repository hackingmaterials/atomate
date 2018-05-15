# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

# This module defines firetasks for writing QChem input files

from atomate.utils.utils import load_class
from fireworks import FiretaskBase, explicit_serialize

__author__ = 'Brandon Wood'
__email__ = "b.wood@berkeley.edu"


@explicit_serialize
class WriteInputFromIOSet(FiretaskBase):
    """
    Writes QChem Input files from input sets. A dictionary is passed to WriteInputFromIOSet where
    parameters are given as keys in the dictionary.

    required_params:
        molecule (Molecule): molecule
        qc_input_set (QChemDictSet or str): Either a QChemDictSet object or a string
        name for the QChem input set (e.g., "OptSet").

    optional_params:
        qchem_input_params (dict): When using a string name for QChem input set, use this as a dict
        to specify kwargs for instantiating the input set parameters. For example, if you want
        to change the DFT_rung, you should provide: {"DFT_rung": ...}.
        This setting is ignored if you provide the full object representation of a QChemDictSet
        rather than a String.
    """

    required_params = ["molecule", "qchem_input_set"]
    optional_params = ["qchem_input_params"]

    def run_task(self, fw_spec):
        # if a full QChemDictSet object was provided
        if hasattr(self['qchem_input_set'], 'write_file'):
            qcin = self['qchem_input_set']

        # if QCInputSet String + parameters was provided
        else:
            qcin_cls = load_class("pymatgen.io.qchem_io.sets", self["qchem_input_set"])
            qcin = qcin_cls(self["molecule"], **self.get("qchem_input_params", {}))
        # we might need to add the filename as a required param
        qcin.write_file("mol.qin")

@explicit_serialize
class WriteInput(FiretaskBase):
    """
    Writes QChem input file from QCInput object.

    required_params:
        qc_input (QCInput): QCInput object

    """
    required_params = ["qc_input"]

    def run_task(self, fw_spec):
        # if a QCInput object is provided
        qcin = self['qc_input']
        qcin.write_file("mol.qin")
