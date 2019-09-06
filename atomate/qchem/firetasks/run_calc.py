# coding: utf-8


# This module defines tasks that support running QChem in various ways.


import shutil
import os
import subprocess

from pymatgen.io.qchem.inputs import QCInput

from custodian import Custodian
from custodian.qchem.handlers import QChemErrorHandler
from custodian.qchem.jobs import QCJob

from fireworks import explicit_serialize, FiretaskBase

from atomate.utils.utils import env_chk, get_logger
import numpy as np

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "5/11/18"
__credits__ = "Shyam Dwaraknath, Xiaohui Qu, Shyue Ping Ong, Anubhav Jain"

logger = get_logger(__name__)


@explicit_serialize
class RunQChemDirect(FiretaskBase):
    """
    Execute a command directly (no custodian).

    Required params:
        qchem_cmd (str): The name of the full command line call to run. This should include any
                         flags for parallelization, saving scratch, and input / output files.
                         Does NOT support env_chk.
    Optional params:
        scratch_dir (str): Path to the scratch directory. Defaults to "/dev/shm/qcscratch/".
                           Supports env_chk.

    """

    required_params = ["qchem_cmd"]
    optional_params = ["scratch_dir"]

    def run_task(self, fw_spec):
        cmd = self.get("qchem_cmd")
        scratch_dir = env_chk(self.get("scratch_dir"), fw_spec)
        if scratch_dir == None:
            scratch_dir = "/dev/shm/qcscratch/"
        os.putenv("QCSCRATCH", scratch_dir)

        logger.info("Running command: {}".format(cmd))
        return_code = subprocess.call(cmd, shell=True)
        logger.info("Command {} finished running with return code: {}".format(
            cmd, return_code))


@explicit_serialize
class RunQChemCustodian(FiretaskBase):
    """
    Run QChem using custodian "on rails", i.e. in a simple way that supports most common options.

    Required params:
        qchem_cmd (str): The name of the full executable for running QChem. Note that this is
                         explicitly different from qchem_cmd in RunQChemDirect because it does
                         not include any flags and should only be the call to the executable.
                         Supports env_chk.

    Optional params:
        multimode (str): Parallelization scheme, either openmp or mpi. Defaults to openmp.
                         Supports env_chk.
        input_file (str): Name of the QChem input file. Defaults to "mol.qin".
        output_file (str): Name of the QChem output file. Defaults to "mol.qout"
        max_cores (int): Maximum number of cores to parallelize over. Supports env_chk.
        qclog_file (str): Name of the file to redirect the standard output to. None means
                          not to record the standard output. Defaults to None.
        suffix (str): String to append to the file in postprocess.
        scratch_dir (str): QCSCRATCH directory. Defaults to "/dev/shm/qcscratch/".
                           Supports env_chk.
        save_scratch (bool): Whether to save scratch directory contents. Defaults to False.
        save_name (str): Name of the saved scratch directory. Defaults to "default_save_name".
        max_errors (int): Maximum # of errors to fix before giving up (default=5)
        job_type (str): Choose from "normal" (default) and "opt_with_frequency_flattener"
        handler_group (str): Group of handlers to use. See handler_groups dict in the code
                             for the groups and complete list of handlers in each group.
        gzip_output (bool): gzip output (default=T)

        *** Just for opt_with_frequency_flattener ***
        max_iterations (int): Number of perturbation -> optimization -> frequency iterations
                              to perform. Defaults to 10.
        max_molecule_perturb_scale (float): The maximum scaled perturbation that can be
                                            applied to the molecule. Defaults to 0.3.

    """
    required_params = ["qchem_cmd"]
    optional_params = [
        "multimode", "input_file", "output_file", "max_cores", "qclog_file",
        "suffix", "scratch_dir", "save_scratch", "save_name", "max_errors",
        "max_iterations", "max_molecule_perturb_scale", "linked",
        "job_type", "handler_group", "gzipped_output"
    ]

    def run_task(self, fw_spec):

        # initialize variables
        qchem_cmd = env_chk(self["qchem_cmd"], fw_spec)
        multimode = env_chk(self.get("multimode"), fw_spec)
        if multimode == None:
            multimode = "openmp"
        input_file = self.get("input_file", "mol.qin")
        output_file = self.get("output_file", "mol.qout")
        max_cores = env_chk(self["max_cores"], fw_spec)
        qclog_file = self.get("qclog_file", "mol.qclog")
        suffix = self.get("suffix", "")
        scratch_dir = env_chk(self.get("scratch_dir"), fw_spec)
        if scratch_dir == None:
            scratch_dir = "/dev/shm/qcscratch/"
        save_scratch = self.get("save_scratch", False)
        save_name = self.get("save_name", "saved_scratch")
        max_errors = self.get("max_errors", 5)
        max_iterations = self.get("max_iterations", 10)
        linked = self.get("linked", False)
        max_molecule_perturb_scale = self.get("max_molecule_perturb_scale",
                                              0.3)
        job_type = self.get("job_type", "normal")
        gzipped_output = self.get("gzipped_output", True)

        handler_groups = {
            "default": [
                QChemErrorHandler(
                    input_file=input_file, output_file=output_file)
            ],
            "no_handler": []
        }

        # construct jobs
        if job_type == "normal":
            jobs = [
                QCJob(
                    qchem_command=qchem_cmd,
                    max_cores=max_cores,
                    multimode=multimode,
                    input_file=input_file,
                    output_file=output_file,
                    qclog_file=qclog_file,
                    suffix=suffix,
                    scratch_dir=scratch_dir,
                    save_scratch=save_scratch,
                    save_name=save_name)
            ]
        elif job_type == "opt_with_frequency_flattener":
            if linked:
                jobs = QCJob.opt_with_frequency_flattener(
                    qchem_command=qchem_cmd,
                    multimode=multimode,
                    input_file=input_file,
                    output_file=output_file,
                    qclog_file=qclog_file,
                    max_iterations=max_iterations,
                    linked=linked,
                    max_cores=max_cores)
            else:
                jobs = QCJob.opt_with_frequency_flattener(
                    qchem_command=qchem_cmd,
                    multimode=multimode,
                    input_file=input_file,
                    output_file=output_file,
                    qclog_file=qclog_file,
                    max_iterations=max_iterations,
                    linked=linked,
                    max_molecule_perturb_scale=max_molecule_perturb_scale,
                    scratch_dir=scratch_dir,
                    save_scratch=save_scratch,
                    save_name=save_name,
                    max_cores=max_cores)

        else:
            raise ValueError("Unsupported job type: {}".format(job_type))

        # construct handlers
        handlers = handler_groups[self.get("handler_group", "default")]

        c = Custodian(
            handlers,
            jobs,
            max_errors=max_errors,
            gzipped_output=gzipped_output)

        c.run()


@explicit_serialize
class RunNoQChem(FiretaskBase):
    """
    Do NOT run QChem. Do nothing.
    """

    def run_task(self, fw_spec):
        pass


@explicit_serialize
class RunQChemFake(FiretaskBase):
    """
     QChem Emulator

     Required params:
         ref_dir (string): Path to reference qchem run directory with input file in the folder
            named "mol.qin" and output file in the folder named "mol.qout".

     """
    required_params = ["ref_dir"]
    optional_params = ["input_file"]

    def run_task(self, fw_spec):
        self._verify_inputs()
        self._clear_inputs()
        self._generate_outputs()

    def _verify_inputs(self):
        input_file = self.get("input_file", "mol.qin")
        user_qin = QCInput.from_file(os.path.join(os.getcwd(), "mol.qin"))

        # Check mol.qin
        ref_qin = QCInput.from_file(os.path.join(self["ref_dir"], input_file))

        np.testing.assert_equal(ref_qin.molecule.species,
                                user_qin.molecule.species)
        np.testing.assert_allclose(
            ref_qin.molecule.cart_coords,
            user_qin.molecule.cart_coords,
            atol=0.0001)
        for key in ref_qin.rem:
            if user_qin.rem.get(key) != ref_qin.rem.get(key):
                raise ValueError("Rem key {} is inconsistent!".format(key))
        if ref_qin.opt is not None:
            for key in ref_qin.opt:
                if user_qin.opt.get(key) != ref_qin.opt.get(key):
                    raise ValueError("Opt key {} is inconsistent!".format(key))
        if ref_qin.pcm is not None:
            for key in ref_qin.pcm:
                if user_qin.pcm.get(key) != ref_qin.pcm.get(key):
                    raise ValueError("PCM key {} is inconsistent!".format(key))
        if ref_qin.solvent is not None:
            for key in ref_qin.solvent:
                if user_qin.solvent.get(key) != ref_qin.solvent.get(key):
                    raise ValueError(
                        "Solvent key {} is inconsistent!".format(key))

        logger.info("RunQChemFake: verified input successfully")

    @staticmethod
    def _clear_inputs():
        p = os.path.join(os.getcwd(), "mol.qin")
        if os.path.exists(p):
            os.remove(p)

    def _generate_outputs(self):
        # pretend to have run QChem by copying pre-generated output from reference dir to cur dir
        for file_name in os.listdir(self["ref_dir"]):
            full_file_name = os.path.join(self["ref_dir"], file_name)
            if os.path.isfile(full_file_name):
                shutil.copy(full_file_name, os.getcwd())
        logger.info("RunQChemFake: ran fake QChem, generated outputs")
