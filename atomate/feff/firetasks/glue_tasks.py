# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

import os

from fireworks import explicit_serialize, FiretaskBase

from atomate.utils.fileio import FileClient
from atomate.common.firetasks.glue_tasks import get_calc_loc, CopyFiles

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

# TODO: @matk86 this has nothing to do with FEFF, apart from certain default files to exclude. It
# also has functionalities that might be useful to other CopyXOutputs style tasks. Can you
# generalize this and put it in a common location, e.g. atomate.common.firetasks.glue_tasks?
# Ideally specific CopyXTasks could subclass the base class with defaults or simply use the base
# class as-is.
# -computron

@explicit_serialize
class CopyFeffOutputs(CopyFiles):
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

        calc_loc = get_calc_loc(self["calc_loc"], fw_spec["calc_locs"]) if self.get("calc_loc") else {}
        exclude_files = self.get("exclude_files", ["feff.inp", "xmu.dat"])

        self.setup(self.get("calc_dir", None), filesystem=self.get("filesystem", None), exclude_files=exclude_files,
                   from_path_dict=calc_loc)
        self.copy()
