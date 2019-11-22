# coding: utf-8


# This module defines tasks that support running QChem in various ways.


import shutil
import os
import subprocess
import logging

from pymatgen.io.qchem.inputs import QCInput
from monty.serialization import loadfn, dumpfn

from custodian import Custodian
from custodian.qchem.handlers import QChemErrorHandler
from custodian.qchem.jobs import QCJob

from monty.tempfile import ScratchDir
from monty.shutil import compress_file, decompress_file
from pymatgen.command_line.critic2_caller import Critic2Output
from fireworks import explicit_serialize, FiretaskBase

from atomate.utils.utils import env_chk#, get_logger
import numpy as np

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "11/20/18"


# logging.basicConfig(filename="critic.log")
# logger = get_logger(__name__)
logger = logging.getLogger(__name__)

@explicit_serialize
class RunCritic2(FiretaskBase):
    """


    """

    required_params = ["molecule", "cube_file"]

    def run_task(self, fw_spec):
        molecule = self.get("molecule")
        cube = self.get("cube_file")

        compress_at_end = False

        if cube[-3:] == ".gz":
            compress_at_end = True
            decompress_file(cube)
            cube = cube[:-3]

        input_script = ["molecule "+cube]
        input_script += ["load "+cube]
        input_script += ["auto"]
        input_script += ["yt"]
        input_script = "\n".join(input_script)

        # with ScratchDir(".") as temp_dir:
        #     os.chdir(temp_dir)
        with open('input_script.cri', 'w') as f:
            f.write(input_script)
        args = ["critic2", "input_script.cri"]

        rs = subprocess.Popen(args,
                              stdout=subprocess.PIPE,
                              stdin=subprocess.PIPE, close_fds=True)

        stdout, stderr = rs.communicate()
        stdout = stdout.decode()

        if stderr:
            stderr = stderr.decode()
            warnings.warn(stderr)

        if rs.returncode != 0:
            # raise RuntimeError("critic2 exited with return code {}.".format(rs.returncode))
            logger.info("critic2 exited with return code {}.".format(rs.returncode))

        logger.info(stderr)
        logger.info(stdout)

        output = Critic2Output(molecule, stdout)
        dumpfn(output.processed_dict,"processed_critic2.json")

        if compress_at_end:
            compress_file(cube)


