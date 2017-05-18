# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines tasks for writing FEFF input sets.
"""

from pymatgen.io.feff.inputs import Paths

from fireworks import FiretaskBase, explicit_serialize

from atomate.utils.utils import load_class

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

# TODO: @matk86 - can you change the format of feff_input_set when a String? Rather than just
# the class name, require the entire path. e.g. "pymatgen.io.feff.sets.MPXANESSet" and put this
# as an example in the parameter documentation for feff_input_set. Thus the code would look like:
# modname, classname = feff_input_set.strip().rsplit(".", 1)
# cls_ = load_class(modname, classname)
# Also update the tests and any workflows, etc.
# This is the convention for other tasks (e.g. ToDbTask) and will avoid confusion. -computron

@explicit_serialize
class WriteFeffFromIOSet(FiretaskBase):
    """
    Generate FEFF input (feff.inp) from the given InputSet object or InputSet name

    Required_params:
        absorbing_atom (str): absorbing atom symbol
        structure (Structure): input structure
        feff_input_set (str or FeffDictSet subclass): The inputset for setting params

    Optional_params:
        radius (float): cluster radius in angstroms
        other_params (dict): **kwargs to pass into the desired InputSet if using str feff_input_set
    """
    required_params = ["absorbing_atom", "structure", "feff_input_set"]
    optional_params = ["radius", "other_params"]

    def run_task(self, fw_spec):
        # if a full FeffInputSet object is provided:
        if hasattr(self['feff_input_set'], 'write_input'):
            fis = self['feff_input_set']

        # else if inputset String + parameters was provided
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
        feff_input_set: (FeffDictSet subclass)
        paths (list): list of paths. A path = list of site indices.

    Optional_params:
        degeneracies (list): list of path degeneracies.
    """
    required_params = ["feff_input_set", "paths"]
    optional_params = ["degeneracies"]

    def run_task(self, fw_spec):
        atoms = self['feff_input_set'].atoms
        paths = Paths(atoms, self["paths"], degeneracies=self.get("degeneracies", []))
        paths.write_file()
