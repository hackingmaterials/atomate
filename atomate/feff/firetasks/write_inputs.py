# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines tasks for writing FEFF input sets.
"""

from fireworks import FiretaskBase, explicit_serialize

from atomate.utils.utils import load_class

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'


@explicit_serialize
class WriteFeffFromIOSet(FiretaskBase):
    """
    Generate FEFF input(feff.inp) from the given inputset object or inputset name

    Required_params:
        absorbing_atom (str): absorbing atom symbol
        structure (Structure): input structure

    Optional_params:
        radius (float): cluster radius in angstroms
        other_params (dict)
    """
    required_params = ["absorbing_atom", "structure", "feff_input_set"]
    optional_params = ["radius", "other_params"]

    def run_task(self, fw_spec):
        # if a full object is provided.
        if hasattr(self['feff_input_set'], 'write_input'):
            fis = self['feff_input_set']

        # if inputset String + parameters was provided
        else:
            fis_cls = load_class("pymatgen.io.feff.sets", self["feff_input_set"])
            fis = fis_cls(self["absorbing_atom"], self["structure"], self.get("radius", 10.0),
                          **self.get("other_params", {}))

        fis.write_input(".")


@explicit_serialize
class WriteEXAFSPaths(FiretaskBase):
    """
    Write the scattering paths to paths.dat file.

    Required_params:
        feff_input_set (FeffDictSet)
        paths (list): list of paths. path = list of site indices.

    Optional_params:
        degeneracies (list): list of path degeneracies.
    """
    required_params = ["feff_input_set", "paths"]
    optional_params = ["degeneracies"]

    def run_task(self, fw_spec):

        from pymatgen.io.feff.inputs import Paths

        atoms = self['feff_input_set'].atoms
        paths = Paths(atoms, self["paths"], degeneracies=self.get("degeneracies", []))
        paths.write_file()
