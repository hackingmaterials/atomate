# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

import os

import six
import re
import monty
import operator

from atomate.utils.utils import env_chk, load_class
from fireworks import explicit_serialize, FiretaskBase, FWAction

__author__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'


@explicit_serialize
class PassCalcLocs(FiretaskBase):
    """
    Passes information about where the current calculation is located
    for the next FireWork. This is achieved by passing a key to
    the fw_spec called "calc_locs" with this information.

    Required params:
        name (str): descriptive name for this calculation file/dir

    Optional params:
        filesystem (str or custom user format): name of filesystem. Supports env_chk. 
            defaults to None
        path (str): The path to the directory containing the calculation. defaults to
            current working directory.
    """

    required_params = ["name"]
    optional_params = ["filesystem", "path"]

    def run_task(self, fw_spec):
        calc_locs = list(fw_spec.get("calc_locs", []))
        calc_locs.append({"name": self["name"],
                          "filesystem": env_chk(self.get('filesystem', None), fw_spec),
                          "path": self.get("path", os.getcwd())})

        return FWAction(mod_spec=[{'_push_all': {'calc_locs': calc_locs}}])


def get_calc_loc(target_name, calc_locs):
    """
    This is a helper method that helps you pick out a certain calculation
    from an array of calc_locs.

    There are three modes:
        - If you set target_name to a String, search for most recent calc_loc
            with matching nameget_
        - Otherwise, return most recent calc_loc overall

    Args:
        target_name: (bool or str) If str, will search for calc_loc with
            matching name, else use most recent calc_loc
        calc_locs: (dict) The dictionary of all calc_locs

    Returns:
        (dict) dict with subkeys path, filesystem, and name
    """

    if isinstance(target_name, six.string_types):
        for doc in reversed(calc_locs):
            if doc["name"] == target_name:
                return doc
        raise ValueError("Could not find the target_name: {}".format(target_name))
    else:
        return calc_locs[-1]


@explicit_serialize
class PassResult(FiretaskBase):
    """
    Passes properties and corresponding user-specified data resulting
    from a run from parent to child fireworks.  Uses a string syntax
    similar to Mongo-style queries to designate values of output
    file dictionaries to retrieve.  For example, one could specify
    a task to pass the stress from the current calculation using:
    
    PassResult(pass_dict={'stress': ">>ionic_steps.-1.stress"})

    Required params:
        pass_dict (dict): dictionary designating keys and values to pass
            to child fireworks.  If value is a string beginning with '>>',
            the firework will search the parsed VASP output dictionary
            for the designated property by following the sequence of keys
            separated with periods, e. g. ">>ionic_steps.-1.stress" is used
            to designate the stress from the last ionic_step. If the value
            is not a string or does not begin with ">>", it is passed as is.
        parse_class (str): string representation of complete path to a class 
            with which to parse the output, e. g. pymatgen.io.vasp.Vasprun
            or pymatgen.io.feff.LDos.from_file, class must be MSONable
        parse_kwargs (str): dict of kwargs for the parse class,
            e. g. {"filename": "vasprun.xml", "parse_dos": False, 
            "parse_eigen": False}

    Optional params:
        calc_dir (str): path to dir that contains VASP output files, defaults
            to '.', e. g. current directory
        mod_spec_cmd (str): command to issue for mod_spec, e. g. "_set" or "_push",
            defaults to "_set"
        mod_spec_key (str): key to pass to mod_spec _set dictmod command, defaults
            to "prev_calc_result"
    """

    required_params = ["pass_dict", "parse_class", "parse_kwargs"]
    optional_params = ["calc_dir", "mod_spec_cmd", "mod_spec_key"]
                       
    def run_task(self, fw_spec):
        pass_dict = self.get("pass_dict")
        parse_kwargs = self.get("parse_kwargs")
        pc_string = self.get("parse_class")
        parse_class = load_class(*pc_string.rsplit(".", 1))
        calc_dir = self.get("calc_dir", ".")
        with monty.os.cd(calc_dir):
            result = parse_class(**parse_kwargs)

        pass_dict = recursive_get_result(pass_dict, result)
        mod_spec_key = self.get("mod_spec_key", "prev_calc_result")
        mod_spec_cmd = self.get("mod_spec_cmd", "_set")
        return FWAction(mod_spec=[{mod_spec_cmd: {mod_spec_key: pass_dict}}])


intfmt = re.compile("[-+]?\d+$")
def recursive_get_result(d, result):
    """
    Helper function for PassResult that gets designated
    keys or values of d (i. e. those that start with "d>>"
    or "a>>") from the corresponding entry in result_dict, 
    similar to FireWorks recursive_deserialize

    Note that the plain ">>" notation will get a key from
    the result.as_dict() object and may use MongoDB
    dot notation, while "a>>" will get an attribute 
    of the object
    """
    if isinstance(d, six.string_types) and d[:2] == ">>":
        # convert integer like strings to ints
        keychain = [int(w) if intfmt.match(w) else w for w in d[2:].split('.')]
        return reduce(operator.getitem, keychain, result.as_dict())

    elif isinstance(d, six.string_types) and d[:3] == "a>>":
        return getattr(result, d[3:])
    
    elif isinstance(d, dict):
        return {k: recursive_get_result(v, result) for k, v in d.items()}
    
    elif isinstance(d, (list, tuple)):
        return [recursive_get_result(i, result) for i in d] 
    
    else:
        return d
