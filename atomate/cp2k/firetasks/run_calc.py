# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

from monty.os.path import zpath
from monty.serialization import loadfn

from atomate.vasp.config import HALF_KPOINTS_FIRST_RELAX

"""
This module defines tasks that support running vasp in various ways.
"""

import shutil
import shlex
import os
import six
import subprocess

from custodian import Custodian
from custodian.cp2k.handlers import UnconvergedScfErrorHandler
from custodian.cp2k.jobs import Cp2kJob
from custodian.vasp.validators import VasprunXMLValidator, VaspFilesValidator

from fireworks import explicit_serialize, FiretaskBase, FWAction

from atomate.utils.utils import env_chk, get_logger
from atomate.vasp.config import CUSTODIAN_MAX_ERRORS

__author__ = 'Nicholas Winner <nwinner@berkeley.edu>'

logger = get_logger(__name__)


@explicit_serialize
class RunCp2KDirect(FiretaskBase):
    """
    Execute a command directly (no custodian).

    Required params:
        cmd (str): the name of the full executable to run. Supports env_chk.
    Optional params:
        expand_vars (str): Set to true to expand variable names in the cmd.
    """

    required_params = ["cp2k_cmd"]
    optional_params = ["expand_vars"]

    def run_task(self, fw_spec):
        cmd = env_chk(self["cp2k_cmd"], fw_spec)
        if self.get("expand_vars", False):
            cmd = os.path.expandvars(cmd)

        logger.info("Running command: {}".format(cmd))
        return_code = subprocess.call(cmd, shell=True)
        logger.info("Command {} finished running with returncode: {}".format(cmd, return_code))


@explicit_serialize
class RunCp2KCustodian(FiretaskBase):
    """
    Run CP2K using custodian "on rails", i.e. in a simple way that supports most common options.

    Required params:
        cp2k_cmd (str): the name of the full executable for running Cp2k. Supports env_chk.

    Optional params:
        job_type: (str) - choose from "normal" (default)
        handler_group: (str or [ErrorHandler]) - group of handlers to use. See handler_groups dict in the code for
            the groups and complete list of handlers in each group. Alternatively, you can
            specify a list of ErrorHandler objects.
        max_force_threshold: (float) - if >0, adds MaxForceErrorHandler. Not recommended for
            nscf runs.
        scratch_dir: (str) - if specified, uses this directory as the root scratch dir.
            Supports env_chk.
        gzip_output: (bool) - gzip output (default=T)
        max_errors: (int) - maximum # of errors to fix before giving up (default=5)
        ediffg: (float) shortcut for setting EDIFFG in special custodian jobs
        auto_npar: (bool) - use auto_npar (default=F). Recommended set to T
            for single-node jobs only. Supports env_chk.
        gamma_vasp_cmd: (str) - cmd for Gamma-optimized VASP compilation.
            Supports env_chk.
        wall_time (int): Total wall time in seconds. Activates WalltimeHandler if set.
    """
    required_params = ["cp2k_cmd"]
    optional_params = []

    def run_task(self, fw_spec):

        handler_groups = {
            "default": [UnconvergedScfErrorHandler()],
            "strict": [],
            "md": [],
            "no_handler": []
            }

        cp2k_cmd = env_chk(self["cp2k_cmd"], fw_spec)

        if isinstance(cp2k_cmd, six.string_types):
            cp2k_cmd = os.path.expandvars(cp2k_cmd)
            cp2k_cmd = shlex.split(cp2k_cmd)

        # initialize variables
        job_type = self.get("job_type", "normal")
        scratch_dir = env_chk(self.get("scratch_dir"), fw_spec)
        gzip_output = self.get("gzip_output", True)
        max_errors = self.get("max_errors", CUSTODIAN_MAX_ERRORS)

        # construct jobs
        if job_type == "normal":
            jobs = [Cp2kJob(cp2k_cmd)]

        # construct handlers
        handler_group = self.get("handler_group", "default")
        if isinstance(handler_group, six.string_types):
            handlers = handler_groups[handler_group]
        else:
            handlers = handler_group

        # construct validators
        validators = []

        c = Custodian(handlers, jobs, validators=validators, max_errors=max_errors,
                      scratch_dir=scratch_dir, gzipped_output=gzip_output)

        c.run()

        if os.path.exists(zpath("custodian.json")):
            stored_custodian_data = {"custodian": loadfn(zpath("custodian.json"))}
            return FWAction(stored_data=stored_custodian_data)
