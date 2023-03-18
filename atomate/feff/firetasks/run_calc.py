"""
This module defines tasks to run FEFF.
"""

import subprocess

from fireworks import FiretaskBase, explicit_serialize

from atomate.utils.utils import env_chk, get_logger

__author__ = "Kiran Mathew"
__email__ = "kmathew@lbl.gov"

logger = get_logger(__name__)


@explicit_serialize
class RunFeffDirect(FiretaskBase):
    """
    Run FEFF.

    Required params:
        feff_cmd (str): the name of the full executable for running FEFF (supports env_chk)
    """

    required_params = ["feff_cmd"]

    def run_task(self, fw_spec):
        feff_cmd = env_chk(self["feff_cmd"], fw_spec)
        logger.info(f"Running FEFF using exe: {feff_cmd}")
        return_code = subprocess.call(feff_cmd, shell=True)
        logger.info(f"FEFF finished running with returncode: {return_code}")
