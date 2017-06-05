# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines firetasks for running lammps
"""

import subprocess

from pymatgen.io.lammps.utils import PackmolRunner
from fireworks import explicit_serialize, FiretaskBase
from atomate.utils.utils import get_logger

__author__ = 'Kiran Mathew'
__email__ = "kmathew@lbl.gov"

logger = get_logger(__name__)


@explicit_serialize
class RunPackmol(FiretaskBase):
    """
    Run packmol.

    Required params:
        molecules (list): list of constituent molecules(Molecule objects)
        packing_config (list): list of dict config settings for each molecule in the
            molecules list. eg: config settings for a single moelcule
            [{"number": 1, "inside box":[0,0,0,100,100,100]}]

    Optional params:
        tolerance (float): packmol tolerance
        filetype (string): input/output structure file type
        control_params: packmol control parameters dictionary. Basically all parameters other
            than structure/atoms
        output_file: output file name. The extension will be adjusted according to the filetype
    """

    required_params = ["molecules", "packing_config"]
    optional_params = ["tolerance", "filetype", "control_params", "output_file"]

    def run_task(self, fw_spec):
        pmr = PackmolRunner(self["molecules"], self["packing_config"],
                            tolerance=self.get("tolerance", 2.0), filetype=self.get("filetype", "xyz"),
                            control_params=self.get("control_params", {"nloop": 1000}),
                            output_file=self.get("output_file", "packed_mol.xyz"))
        logger.info("Running packmol")
        pmr.run()
        logger.info("Packmol finished running.")


@explicit_serialize
class RunLammpsDirect(FiretaskBase):
    """
    Run LAMMPS directly (no custodian).

    Required params:
        lammsps_cmd (str): full command string
    """

    required_params = ["lammps_cmd"]

    def run_task(self, fw_spec):
        lammps_cmd = self["lammps_cmd"]
        logger.info("Running LAMMPS using exe: {}".format(lammps_cmd))
        return_code = subprocess.call(lammps_cmd, shell=True)
        logger.info("LAMMPS finished running with returncode: {}".format(return_code))
