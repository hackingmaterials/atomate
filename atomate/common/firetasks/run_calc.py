from custodian import Custodian
from fireworks import explicit_serialize, FiretaskBase


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