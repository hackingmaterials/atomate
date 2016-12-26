# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

import os

from fireworks import explicit_serialize, FiretaskBase, FWAction
from atomate.utils.utils import env_chk

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


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