import subprocess
from fireworks import explicit_serialize, FireTaskBase

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

@explicit_serialize
class RunVaspDirect(FireTaskBase):
    """
    Run VASP directly (no custodian).

    Required params:
        (none)

    Optional params:
        vasp_exe (str): the name of the full executable for running VASP. Used if vasp_exe_fwenv_key not set.
        vasp_exe_env_key (str): name of key in fw_env that contains the name of the full VASP executable.
                                Used to tailor vasp exe for different workers.
    """

    optional_params = ["vasp_exe", "vasp_exe_env_key"]

    def run_task(self, fw_spec):
        # Establish VASP executable
        if "vasp_exe_env_key" in self:
            vasp_exe = fw_spec["_fw_env"][self["vasp_exe_env_key"]]
        elif "vasp_exe" in self:
            vasp_exe = self["vasp_exe"]
        else:
            vasp_exe = "vasp"

        # Run VASP
        print("Running VASP using exe: {}".format(vasp_exe))
        return_code = subprocess.call(vasp_exe, shell=True)
        print("VASP finished running with returncode: {}".format(return_code))


@explicit_serialize
class RunVaspCustodian(FireTaskBase):
    """
    Run VASP using custodian.

    Required params:
        (none)

    Optional params:
        vasp_exe (str): the name of the full executable for running VASP. Used if vasp_exe_fwenv_key not set.
        vasp_exe_env_key (str): name of key in fw_env that contains the name of the full VASP executable.
                                Used to tailor vasp exe for different workers.
    """

    optional_params = ["vasp_exe", "vasp_exe_env_key"]

    def run_task(self, fw_spec):
        # Establish VASP executable
        if "vasp_exe_env_key" in self:
            vasp_exe = fw_spec["_fw_env"][self["vasp_exe_env_key"]]
        elif "vasp_exe" in self:
            vasp_exe = self["vasp_exe"]
        else:
            vasp_exe = "vasp"

        # Run VASP
        pass
