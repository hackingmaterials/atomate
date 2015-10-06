import subprocess
from custodian import Custodian
from fireworks import explicit_serialize, FireTaskBase
from matmethods.utils.utils import env_chk

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

@explicit_serialize
class RunVaspDirect(FireTaskBase):
    """
    Run VASP directly (no custodian).

    Required params:
        vasp_exe (str): the name of the full executable for running VASP. Supports env_chk.
    """

    required_params = ["vasp_exe"]

    def run_task(self, fw_spec):
        vasp_exe = env_chk(self["vasp_exe"], fw_spec)

        print("Running VASP using exe: {}".format(vasp_exe))
        return_code = subprocess.call(vasp_exe, shell=True)
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
        jobs: ([Job]) - a list of custodian jobs to run
        handlers: ([ErrorHandler]) - a list of error handlers

    Optional params:
        job_type: (str) - choose from "normal" (default), "double_relaxation", and "iterative_relaxation"
        scratch_dir: (str) - if specified, uses this directory as the root scratch dir. Supports env_chk.

    """

    required_params = ["jobs", "handlers"]
    optional_params = ["validators", "custodian_params"]

    def run_task(self, fw_spec):
        c = Custodian(self["handlers"], self["jobs"], self.get("validators"), **self.get("custodian_params", {}))
        output = c.run()