# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

import os

import six
import re
import monty
import operator
from functools import reduce

from atomate.utils.utils import env_chk, load_class, recursive_get_result
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


@explicit_serialize
class GrabFilesFromCalcLoc(FiretaskBase):
    """
    Based on CopyVaspOutputs but for general file copying

    This FireTask can be inhertited to use the get_file functionality.

    May need to have arguements to get_file that have the calc_dir of
    the desired calculation file, so this can be generically used by
    whatever class inherits it. This way, an inheriting class can
    decide which variables need to be pulled from fw_spec to determine
    which files need to be copied.
    """

    required_params = ["filenames", "name_prepend","name_append"]
    optional_params = ["calc_dir", "calc_loc"]

    def run_task(self,fw_spec=None):

        if self.get('calc_dir',False):
            filesystem = None
        elif self.get('calc_loc',False):
            calc_loc = get_calc_loc(self.get('calc_loc',False), fw_spec["calc_locs"])
            calc_dir = calc_loc["path"]
            filesystem = calc_loc["filesystem"]
        else:
            raise ValueError("Must specify either calc_dir or calc_loc!")

        fileclient = FileClient(filesystem=filesystem)
        calc_dir = fileclient.abspath(calc_dir)

        all_files = fileclient.listdir(calc_dir)

        if self.get('filenames',False):
            if type(self.get('filenames',False)) == list:
                files_to_copy = self.get('filenames',False)
            elif isinstance(self.get('filenames',False), six.string_types):
                files_to_copy = [ self.get('filenames',False) ]
            else:
                ValueError("Must have a list of strings or a strings!")
        else:
            raise ValueError("Must have a list of filenames!")

        # determine what files need to be copied
        if "$ALL" in self.get('filenames',False):
            files_to_copy = all_files

        # start file copy
        for f in files_to_copy:
            prev_path_full = os.path.join(calc_dir, f)
            # prev_path = os.path.join(os.path.split(calc_dir)[1], f)
            dest_fname = self.get('name_prepend',"")+f+self.get('name_append',"")
            dest_path = os.path.join(os.getcwd(), dest_fname)

            # copy the file (minus the relaxation extension)
            fileclient.copy(prev_path_full, dest_path)
            print(prev_path_full)
            print(dest_path)


@explicit_serialize
class CreateFolder(FiretaskBase):
    """
    FireTask to create new folder with the option of changing directory to the new folder.
    """
    required_params = ["folder_name","change_to"]

    def run_task(self, fw_spec):
        folder = "/"+self.get("folder_name", "test")
        if not os.path.exists(os.getcwd() + folder):
            os.makedirs(os.getcwd() + folder)
        if self.get("change_to", False):
            os.chdir(os.getcwd() + folder)

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
