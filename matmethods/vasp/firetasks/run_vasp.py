import shlex
import subprocess
from custodian import Custodian
from custodian.vasp.jobs import VaspJob
from fireworks import explicit_serialize, FireTaskBase
from matmethods.utils.utils import env_chk

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

@explicit_serialize
class RunVaspDirect(FireTaskBase):
    """
    Run VASP directly (no custodian).

    Required params:
        vasp_cmd (str): the name of the full executable for running VASP. Supports env_chk.
    """

    required_params = ["vasp_cmd"]

    def run_task(self, fw_spec):
        vasp_cmd = env_chk(self["vasp_cmd"], fw_spec)

        print("Running VASP using exe: {}".format(vasp_cmd))
        return_code = subprocess.call(vasp_cmd, shell=True)
        print("VASP finished running with returncode: {}".format(return_code))


@explicit_serialize
class RunVaspCustodianFromObjects(FireTaskBase):
    """
    Run VASP using custodian in a generic manner using built-in custodian objects

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
        c = Custodian(self["handlers"], self["jobs"], self.get("validators"), **self.get("custodian_params", {}))
        output = c.run()

@explicit_serialize
class RunVaspCustodianOnRails(FireTaskBase):
    """
    Run VASP using custodian in a generic manner using built-in custodian objects

    Required params:
        vasp_cmd (str): the name of the full executable for running VASP. Supports env_chk.

    Optional params:
        job_type: (str) - choose from "normal" (default), "double_relaxation", and "iterative_relaxation"
        scratch_dir: (str) - if specified, uses this directory as the root scratch dir. Supports env_chk.
        gzip_output: (bool) - gzip output (default=T)
        max_errors: (int) - max errors (default=2)
        auto_npar: (bool) - use auto_npar (default=F). Supports env_chk.
        gamma_vasp_cmd: (str) - cmd for Gamma-optimized VASP compilation. Supports env_chk.

    """

    required_params = ["vasp_cmd"]
    optional_params = ["job_type", "scratch_dir", "gzip_output", "max_errors", "auto_npar", "gamma_vasp_cmd"]

    def run_task(self, fw_spec):
        vasp_cmd = env_chk(self["vasp_cmd"], fw_spec)
        if isinstance(vasp_cmd, basestring):
            vasp_cmd = shlex.split(vasp_cmd)

        # initialize variables

        job_type = self.get("job_type", "normal")
        scratch_dir = env_chk(self.get("scratch_dir"))
        gzip_output = self.get("gzip_output", True)
        max_errors = self.get("max_errors", 2)
        auto_npar = self.get("auto_npar", True)
        gamma_vasp_cmd = self.get("gamma_vasp_cmd")

        # construct jobs
        jobs = []
        if job_type == "normal":
            jobs = [VaspJob(vasp_cmd, gzipped=gzip_output, default_vasp_input_set=None, auto_npar=auto_npar, gamma_vasp_cmd=gamma_vasp_cmd)]
        elif job_type == "double_relaxation":
            jobs = VaspJob.double_relaxation_run(vasp_cmd, gzipped=gzip_output)
        elif job_type == "iterative_relaxation":
            jobs = VaspJob.full_opt_run(vasp_cmd, auto_npar=auto_npar, max_steps=5)
        else:
            raise ValueError("Unsupported job type: {}".format(job_type))


        c = Custodian(self["handlers"], jobs, self.get("validators"), **self.get("custodian_params", {}))
        output = c.run()