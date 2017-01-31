# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines tasks that support running vasp in various ways.
"""

import shutil
import shlex
import subprocess
import os
import six

from pymatgen.io.vasp import Incar, Kpoints, Poscar, Potcar
from pymatgen.io.vasp.sets import get_vasprun_outcar
from pymatgen.electronic_structure.boltztrap import BoltztrapRunner

from custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler, AliasingErrorHandler, MeshSymmetryErrorHandler, \
    UnconvergedErrorHandler, MaxForceErrorHandler, PotimErrorHandler, FrozenJobErrorHandler, \
    NonConvergingErrorHandler, PositiveEnergyErrorHandler, WalltimeHandler
from custodian.vasp.jobs import VaspJob
from custodian.vasp.validators import VasprunXMLValidator, VaspFilesValidator

from fireworks import explicit_serialize, FiretaskBase

from atomate.utils.utils import env_chk, get_logger

__author__ = 'Anubhav Jain <ajain@lbl.gov>'
__credits__ = 'Shyue Ping Ong <ong.sp>'

logger = get_logger(__name__)


@explicit_serialize
class RunVaspDirect(FiretaskBase):
    """
    Run VASP directly (no custodian).

    Required params:
        vasp_cmd (str): the name of the full executable for running VASP.
        Supports env_chk.
    """

    required_params = ["vasp_cmd"]

    def run_task(self, fw_spec):
        vasp_cmd = env_chk(self["vasp_cmd"], fw_spec)
        logger.info("Running VASP using exe: {}".format(vasp_cmd))
        return_code = subprocess.call(vasp_cmd, shell=True)
        logger.info("VASP finished running with returncode: {}".format(return_code))


@explicit_serialize
class RunVaspCustodianFromObjects(FiretaskBase):
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


@explicit_serialize
class RunVaspCustodian(FiretaskBase):
    """
    Run VASP using custodian "on rails", i.e. in a simple way that supports
    most common options.

    Required params:
        vasp_cmd (str): the name of the full executable for running VASP.
            Supports env_chk.

    Optional params:
        job_type: (str) - choose from "normal" (default),
            "double_relaxation_run" (two consecutive jobs), and "full_opt_run"
        handler_group: (str) - group of handlers to use. See handler_groups
            dict in the code for the groups and complete list of handlers in
            each group.
        max_force_threshold: (float) - if >0, adds MaxForceErrorHandler.
            Not recommended for nscf runs.
        scratch_dir: (str) - if specified, uses this directory as the root
            scratch dir. Supports env_chk.
        gzip_output: (bool) - gzip output (default=T)
        max_errors: (int) - maximum # of errors to fix before giving up (default=5)
        ediffg: (float) shortcut for setting EDIFFG in special custodian jobs
        auto_npar: (bool) - use auto_npar (default=F). Recommended set to T
            for single-node jobs only. Supports env_chk.
        gamma_vasp_cmd: (str) - cmd for Gamma-optimized VASP compilation.
            Supports env_chk.
        wall_time (int): Total wall time in seconds. Activates WalltimeHandler if set.
    """
    required_params = ["vasp_cmd"]
    optional_params = ["job_type", "handler_group", "max_force_threshold", "scratch_dir",
                       "gzip_output", "max_errors", "ediffg", "auto_npar", "gamma_vasp_cmd",
                       "wall_time"]

    def run_task(self, fw_spec):

        handler_groups = {
            "default": [VaspErrorHandler(), MeshSymmetryErrorHandler(), UnconvergedErrorHandler(),
                        NonConvergingErrorHandler(),PotimErrorHandler(), PositiveEnergyErrorHandler(),
                        FrozenJobErrorHandler()],
            "strict": [VaspErrorHandler(), MeshSymmetryErrorHandler(), UnconvergedErrorHandler(),
                       NonConvergingErrorHandler(),PotimErrorHandler(), PositiveEnergyErrorHandler(),
                       FrozenJobErrorHandler(), AliasingErrorHandler()],
            "md": [VaspErrorHandler(), NonConvergingErrorHandler()],
            "no_handler": []
            }

        vasp_cmd = env_chk(self["vasp_cmd"], fw_spec)
        if isinstance(vasp_cmd, six.string_types):
            vasp_cmd = os.path.expandvars(vasp_cmd)
            vasp_cmd = shlex.split(vasp_cmd)

        # initialize variables
        job_type = self.get("job_type", "normal")
        scratch_dir = env_chk(self.get("scratch_dir"), fw_spec)
        gzip_output = self.get("gzip_output", True)
        max_errors = self.get("max_errors", 5)
        auto_npar = env_chk(self.get("auto_npar"), fw_spec, strict=False, default=False)
        gamma_vasp_cmd = env_chk(self.get("gamma_vasp_cmd"), fw_spec, strict=False, default=None)
        if gamma_vasp_cmd:
            gamma_vasp_cmd = shlex.split(gamma_vasp_cmd)

        # construct jobs
        if job_type == "normal":
            jobs = [VaspJob(vasp_cmd, auto_npar=auto_npar, gamma_vasp_cmd=gamma_vasp_cmd)]
        elif job_type == "double_relaxation_run":
            jobs = VaspJob.double_relaxation_run(vasp_cmd, auto_npar=auto_npar, ediffg=self.get("ediffg"),
                                                 half_kpts_first_relax=False)
        elif job_type == "full_opt_run":
            jobs = VaspJob.full_opt_run(vasp_cmd, auto_npar=auto_npar, ediffg=self.get("ediffg"),
                                        max_steps=5, half_kpts_first_relax=False)
        else:
            raise ValueError("Unsupported job type: {}".format(job_type))

        # construct handlers
        handlers = handler_groups[self.get("handler_group", "default")]

        if self.get("max_force_threshold"):
            handlers.append(MaxForceErrorHandler(max_force_threshold=self["max_force_threshold"]))

        if self.get("wall_time"):
            handlers.append(WalltimeHandler(wall_time=self["wall_time"]))

        validators = [VasprunXMLValidator(), VaspFilesValidator()]

        c = Custodian(handlers, jobs, validators=validators, max_errors=max_errors,
                      scratch_dir=scratch_dir, gzipped_output=gzip_output)

        c.run()


@explicit_serialize
class RunBoltztrap(FiretaskBase):
    """
    Run Boltztrap directly. Requires vasprun.xml and OUTCAR to be
    in current dir.

    Required params:
        (none)

    Optional params:
        scissor: (float) scissor band gap amount in eV (i.e. new gap == scissor)
        tmax: (float) max temperature to evaluate (default = 1300K)
        tgrid: (float) temperature interval (default = 50K)
        doping: ([float]) doping levels you want to compute
        soc: (bool) whether the band structure is calculated with spin-orbit coupling or not
    """

    optional_params = ["scissor", "tmax", "tgrid", "doping", "soc"]

    def run_task(self, fw_spec):
        scissor = self.get("scissor", 0.0)
        tmax = self.get("tmax", 1300)
        tgrid = self.get("tgrid", 50)
        doping = self.get("doping", None)
        soc = self.get("soc", False)

        vasprun, outcar = get_vasprun_outcar(".", parse_dos=True, parse_eigen=True)
        bs = vasprun.get_band_structure()
        nelect = outcar.nelect
        runner = BoltztrapRunner(bs, nelect, scissor=scissor, doping=doping, tmax=tmax, tgrid=tgrid, soc=soc)
        runner.run(path_dir=os.getcwd())


@explicit_serialize
class RunVaspFake(FiretaskBase):
    """
     Vasp Emulator

     Required params:
         ref_dir (string): Path to reference vasp run directory with input files
            in the folder named 'inputs' and output files in the folder named 'outputs'.

     Optional params:
         params_to_check (list): optional list of incar parameters to check.
     """
    required_params = ["ref_dir"]
    optional_params = ["params_to_check"]

    def run_task(self, fw_spec):
        self._verify_inputs()
        self._clear_inputs()
        self._generate_outputs()

    def _verify_inputs(self):
        user_incar = Incar.from_file(os.path.join(os.getcwd(), "INCAR"))
        ref_incar = Incar.from_file(os.path.join(self["ref_dir"], "inputs", "INCAR"))

        # perform some BASIC tests

        # check INCAR
        params_to_check = self.get("params_to_check", [])
        defaults = {"ISPIN": 1, "ISMEAR": 1, "SIGMA": 0.2}
        for p in params_to_check:
            if user_incar.get(p, defaults.get(p)) != ref_incar.get(p, defaults.get(p)):
                raise ValueError("INCAR value of {} is inconsistent!".format(p))

        # check KPOINTS
        user_kpoints = Kpoints.from_file(os.path.join(os.getcwd(), "KPOINTS"))
        ref_kpoints = Kpoints.from_file(os.path.join(self["ref_dir"], "inputs", "KPOINTS"))
        if user_kpoints.style != ref_kpoints.style or user_kpoints.num_kpts != ref_kpoints.num_kpts:
            raise ValueError("KPOINT files are inconsistent! Paths are:\n{}\n{}".format(
                os.getcwd(), os.path.join(self["ref_dir"], "inputs")))

        # check POSCAR
        user_poscar = Poscar.from_file(os.path.join(os.getcwd(), "POSCAR"))
        ref_poscar = Poscar.from_file(os.path.join(self["ref_dir"], "inputs", "POSCAR"))
        if user_poscar.natoms != ref_poscar.natoms or user_poscar.site_symbols != ref_poscar.site_symbols:
            raise ValueError("POSCAR files are inconsistent! Paths are:\n{}\n{}".format(
                os.getcwd(), os.path.join(self["ref_dir"], "inputs")))

        # check POTCAR
        user_potcar = Potcar.from_file(os.path.join(os.getcwd(), "POTCAR"))
        ref_potcar = Potcar.from_file(os.path.join(self["ref_dir"], "inputs", "POTCAR"))
        if user_potcar.symbols != ref_potcar.symbols:
            raise ValueError("POTCAR files are inconsistent! Paths are:\n{}\n{}".format(
                os.getcwd(), os.path.join(self["ref_dir"], "inputs")))
        logger.info("RunVaspFake: verified inputs successfully")

    def _clear_inputs(self):
        for x in ["INCAR", "KPOINTS", "POSCAR", "POTCAR", "CHGCAR", "OUTCAR", "vasprun.xml"]:
            p = os.path.join(os.getcwd(), x)
            if os.path.exists(p):
                os.remove(p)

    def _generate_outputs(self):
        output_dir = os.path.join(self["ref_dir"], "outputs")
        for file_name in os.listdir(output_dir):
            full_file_name = os.path.join(output_dir, file_name)
            if os.path.isfile(full_file_name):
                shutil.copy(full_file_name, os.getcwd())
        logger.info("RunVaspFake: ran fake VASP, generated outputs")
