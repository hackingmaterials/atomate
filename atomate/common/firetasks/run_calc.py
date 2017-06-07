# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

import os
import subprocess

from custodian import Custodian

from fireworks import explicit_serialize, FiretaskBase

from atomate.utils.utils import env_chk, get_logger

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

logger = get_logger(__name__)


@explicit_serialize
class RunCommand(FiretaskBase):
    """
    Execute a command directly (no custodian).

    Required params:
        cmd (str): the name of the full executable to run. Supports env_chk.
    Optional params:
        expand_vars (str): Set to true to expand variable names in the cmd.
    """

    required_params = ["cmd"]
    optional_params = ["expand_vars"]

    def run_task(self, fw_spec):
        cmd = env_chk(self["cmd"], fw_spec)
        if self.get("expand_vars", False):
            cmd = os.path.expandvars(cmd)

        logger.info("Running command: {}".format(cmd))
        return_code = subprocess.call(cmd, shell=True)
        logger.info("Command {} finished running with returncode: {}".format(cmd, return_code))


@explicit_serialize
class RunCustodianFromObjects(FiretaskBase):
    """
    Run VASP using custodian in a generic manner using built-in custodian
    objects

    Required params:
        jobs: ([Job]) - a list of custodian jobs to run
        handlers: ([ErrorHandler]) - a list of error handlers

    Optional params:
        validators: ([Validator]) - a list of Validators
        custodian_params ({}) - dict of all other custodian parameters
    """

    required_params = ["jobs", "handlers"]
    optional_params = ["validators", "custodian_params"]

    def run_task(self, fw_spec):
        c = Custodian(self["handlers"], self["jobs"], self.get("validators"),
                      **self.get("custodian_params", {}))
        c.run()
