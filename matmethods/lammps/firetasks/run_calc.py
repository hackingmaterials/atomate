# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines firetasks for running lammps
"""

import subprocess

from fireworks import explicit_serialize, FireTaskBase


__author__ = 'Kiran Mathew'


@explicit_serialize
class RunLammpsDirect(FireTaskBase):
    """
    Run LAMMPS directly (no custodian).

    Required params:
        lammsps_cmd (str): full command string
    """

    required_params = ["lammps_cmd"]

    def run_task(self, fw_spec):
        lammps_cmd = self["lammps_cmd"]
        print("Running LAMMPS using exe: {}".format(lammps_cmd))
        return_code = subprocess.call(lammps_cmd, shell=True)
        print("LAMMPS finished running with returncode: {}".format(return_code))


@explicit_serialize
class RunLammpsCustodian(FireTaskBase):
    pass