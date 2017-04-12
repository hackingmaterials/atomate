# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

import os

from fireworks import explicit_serialize, FiretaskBase

from atomate.utils.fileio import FileClient
from atomate.common.firetasks.glue_tasks import get_calc_loc

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

# TODO: @matk86 this has nothing to do with FEFF, apart from certain default files to exclude. It
# also has functionalities that might be useful to other CopyXOutputs style tasks. Can you
# generalize this and put it in a common location, e.g. atomate.common.firetasks.glue_tasks?
# Ideally specific CopyXTasks could subclass the base class with defaults or simply use the base
# class as-is.
# -computron

@explicit_serialize
class CopyFeffOutputs(FiretaskBase):
    """
    Copy files from a previous run directory to the current directory.
    Note: must specify either "calc_loc" or "calc_dir" to indicate the directory
        containing the files to copy.

    Optional params:
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
            search for the most recent calc_loc with the matching name.
        calc_dir (str): path to dir that contains VASP output files.
        filesystem (str): remote filesystem. e.g. username@host
        exclude_files (list): list fo filenames to be excluded when copying.
    """

    optional_params = ["calc_loc", "calc_dir", "filesystem", "exclude_files"]

    def run_task(self, fw_spec):

        exclude_files = self.get("exclude_files", ["feff.inp", "xmu.dat"])

        if self.get("calc_dir"):
            calc_dir = self["calc_dir"]
            filesystem = None
        elif self.get("calc_loc"):
            calc_loc = get_calc_loc(self["calc_loc"], fw_spec["calc_locs"])
            calc_dir = calc_loc["path"]
            filesystem = calc_loc["filesystem"]
        else:
            raise ValueError("Must specify either calc_dir or calc_loc!")

        fileclient = FileClient(filesystem=filesystem)
        calc_dir = fileclient.abspath(calc_dir)

        files_to_copy = [f for f in fileclient.listdir(calc_dir) if f not in exclude_files]

        for f in files_to_copy:
            prev_path_full = os.path.join(calc_dir, f)
            dest_path = os.path.join(os.getcwd(), f)
            fileclient.copy(prev_path_full, dest_path)
