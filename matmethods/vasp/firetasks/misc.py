import os

from fireworks import explicit_serialize, FireTaskBase, FWAction
from matmethods.utils.utils import env_chk

__author__ = 'Anubhav Jain <ajain@lbl.gov>'



@explicit_serialize
class PassVaspDirs(FireTaskBase):
    """
    Passes the vasp_locs key. Should be called in the same FireWork as a VASP run.
    Needed for certain downstream FireTasks

    Required params:
        name: descriptive name for this VASP file/dir

    Optional params:
        filesystem: name of filesystem. Supports env_chk. defaults to None
        path: The path to the directory containing the VASP run. defaults to current working directory.
    """

    required_params = ["name"]
    optional_params = ["filesystem", "vaspdir"]

    def run_task(self, fw_spec):
        doc = {"name": self["name"],
               "filesystem": env_chk(self.get('filesystem', None), fw_spec),
               "path": self.get("path", os.getcwd())}

        return FWAction(mod_spec=[{'_push': {'vasp_locs': doc}}])