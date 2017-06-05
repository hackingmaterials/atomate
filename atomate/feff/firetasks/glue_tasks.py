# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

from fireworks import explicit_serialize

from atomate.common.firetasks.glue_tasks import get_calc_loc, CopyFiles

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'


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

        self.setup_copy(self.get("calc_dir", None), filesystem=self.get("filesystem", None),
                        exclude_files=exclude_files, from_path_dict=calc_loc)
        self.copy_files()
