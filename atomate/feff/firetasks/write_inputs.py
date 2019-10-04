# coding: utf-8

"""
This module defines tasks for writing FEFF input sets.
"""

from pymatgen.io.feff.inputs import Paths

from fireworks import FiretaskBase, explicit_serialize

from atomate.utils.utils import load_class

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'


@explicit_serialize
class WriteFeffFromIOSet(FiretaskBase):
    """
    Generate FEFF input (feff.inp) from the given InputSet object or InputSet name

    Required_params:
        absorbing_atom (str): absorbing atom symbol
        structure (Structure): input structure
        feff_input_set (str or FeffDictSet subclass): The inputset for setting params. If string
            then either the entire path to the class or the spectrum type must be provided
            e.g. "pymatgen.io.feff.sets.MPXANESSet" or "XANES"

    Optional_params:
        radius (float): cluster radius in angstroms
        other_params (dict): **kwargs to pass into the desired InputSet if using str feff_input_set
    """
    required_params = ["absorbing_atom", "structure", "feff_input_set"]
    optional_params = ["radius", "other_params"]

    def run_task(self, fw_spec):
        feff_input_set = get_feff_input_set_obj(self["feff_input_set"], self["absorbing_atom"],
                                                self["structure"], self.get("radius", 10.0),
                                                **self.get("other_params", {}))
        feff_input_set.write_input(".")


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


def get_feff_input_set_obj(fis, *args, **kwargs):
    """
    returns feff input set object.

    Args:
        fis (str or FeffDictSet subclass): The inputset for setting params. If string then
            the entire path to the class or the spectrum type must be provided
            e.g. "pymatgen.io.feff.sets.MPXANESSet" or "XANES"
        args (tuple): feff input set args
        kwargs (dict): feff input set kwargs

    Returns:
        FeffDictSet object
    """
    # e.g. "pymatgen.io.feff.sets.MPXANESSet" or "XANES"
    if isinstance(fis, str):
        fis_ = "pymatgen.io.feff.sets.MP{}Set".format(fis) if "pymatgen" not in fis else fis
        modname, classname = fis_.strip().rsplit(".", 1)
        fis_cls = load_class(modname, classname)
        return fis_cls(*args, **kwargs)
    else:
        return fis