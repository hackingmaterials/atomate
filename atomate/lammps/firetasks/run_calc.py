# coding: utf-8


"""
This module defines firetasks for running lammps
"""

import os
import shutil

# from pymatgen.io.lammps.utils import PackmolRunner, LammpsRunner

from fireworks import explicit_serialize, FiretaskBase, FWAction

from atomate.utils.utils import get_logger

__author__ = 'Kiran Mathew'
__email__ = "kmathew@lbl.gov"

logger = get_logger(__name__)


@explicit_serialize
class RunLammpsDirect(FiretaskBase):
    """
    Run LAMMPS directly (no custodian).

    Required params:
        lammsps_cmd (str): lammps command to run sans the input file name.
            e.g. 'mpirun -n 4 lmp_mpi'
    """

    required_params = ["lammps_cmd", "input_filename"]

    def run_task(self, fw_spec):
        lammps_cmd = self["lammps_cmd"]
        input_filename = self["input_filename"]
        lmps_runner = LammpsRunner(input_filename, lammps_cmd)
        stdout, stderr = lmps_runner.run()
        logger.info("LAMMPS finished running: {} \n {}".format(stdout, stderr))


@explicit_serialize
class RunLammpsFake(FiretaskBase):
    """
    Pretend run i.e just copy files from existing run dir.

    Required params:
        ref_dir (str): path to the reference dir
    """

    required_params = ["ref_dir"]

    def run_task(self, fw_spec):
        output_dir = os.path.abspath(self["ref_dir"])
        for file_name in os.listdir(output_dir):
            full_file_name = os.path.join(output_dir, file_name)
            if os.path.isfile(full_file_name):
                shutil.copy(full_file_name, os.getcwd())
        logger.info("Ran fake LAMMPS.")


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
        control_params (dict): packmol control parameters dictionary. Basically all parameters other
            than structure/atoms
        output_file (str): output file name. The extension will be adjusted according to the filetype
        copy_to_current_on_exit (bool): whether or not to copy the packed molecule output file to
            the current directory.
        site_property (str): the specified site property will be restored for the final Molecule object.
    """

    required_params = ["molecules", "packing_config", "packmol_cmd"]
    optional_params = ["tolerance", "filetype", "control_params", "output_file",
                       "copy_to_current_on_exit", "site_property"]

    def run_task(self, fw_spec):
        pmr = PackmolRunner(self["molecules"], self["packing_config"],
                            tolerance=self.get("tolerance", 2.0),
                            filetype=self.get("filetype", "xyz"),
                            control_params=self.get("control_params", {"nloop": 1000}),
                            output_file=self.get("output_file", "packed_mol.xyz"),
                            bin=self["packmol_cmd"])
        logger.info("Running {}".format(self["packmol_cmd"]))
        packed_mol = pmr.run(self.get("copy_to_current_on_exit", False), site_property=self.get("site_property", None))
        logger.info("Packmol finished running.")
        return FWAction(mod_spec=[{'_set': {'packed_mol': packed_mol}}])
