# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

import os

from atomate.utils.utils import env_chk
from atomate.utils.fileio import FileClient
from fireworks import explicit_serialize, FiretaskBase, FWAction
from atomate.utils.utils import get_calc_loc

__author__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'


@explicit_serialize
class PassCalcLocs(FiretaskBase):
    """
    Passes the calc_locs key. Should be called in the same FireWork as a
    the calculation. This passes information about where the current run is located
    for the next FireWork.

    Required params:
        name: descriptive name for this calculation file/dir

    Optional params:
        filesystem: name of filesystem. Supports env_chk. defaults to None
        path: The path to the directory containing the calculation. defaults to
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

    May need to have arguements to get_file that have the calc_dir of the desired calculation file, so this can be
    generically used by whatever class inherits it. This way, an inheriting class can decide which variables need to be
    pulled from fw_spec to determine which files need to be copied.
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