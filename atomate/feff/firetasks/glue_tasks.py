# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

import os

from fireworks import explicit_serialize, FiretaskBase

from atomate.utils.fileio import FileClient
from atomate.utils.utils import get_calc_loc

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'


@explicit_serialize
class CopyFeffOutputs(FiretaskBase):
    """
    Copy files from a previous run directory to the current directory.
    Note: must specify either "calc_loc" or "calc_dir" to indicate the
        directory containing the files to copy.

    Required params:
        (none) - but you must specify either "calc_loc" OR "calc_dir"

    Optional params:
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
            search for the most recent calc_loc with the matching name
        calc_dir (str): path to dir that contains VASP output files.
        filesystem (str): remote filesystem. e.g. username@host
    """

    optional_params = ["calc_loc", "calc_dir", "filesystem"]

    def run_task(self, fw_spec):

        if self.get("calc_dir"):  # direct setting of calc dir - no calc_locs or filesystem!
            calc_dir = self["calc_dir"]
            filesystem = None
        elif self.get("calc_loc"):  # search for calc dir and filesystem within calc_locs
            calc_loc = get_calc_loc(self["calc_loc"], fw_spec["calc_locs"])
            calc_dir = calc_loc["path"]
            filesystem = calc_loc["filesystem"]
        else:
            raise ValueError("Must specify either calc_dir or calc_loc!")

        fileclient = FileClient(filesystem=filesystem)
        calc_dir = fileclient.abspath(calc_dir)

        files_to_copy = [f for f in fileclient.listdir(calc_dir) if f not in ["feff.inp", "xmu.dat"]]

        # start file copy
        for f in files_to_copy:
            prev_path_full = os.path.join(calc_dir, f)
            # prev_path = os.path.join(os.path.split(calc_dir)[1], f)
            dest_path = os.path.join(os.getcwd(), f)

            # copy the file (minus the relaxation extension)
            fileclient.copy(prev_path_full, dest_path)