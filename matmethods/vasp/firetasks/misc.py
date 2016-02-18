import os
import shutil

from fireworks import explicit_serialize, FireTaskBase, FWAction
from matmethods.utils.utils import env_chk

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


@explicit_serialize
class PassVaspLocs(FireTaskBase):
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
    optional_params = ["filesystem", "path"]

    def run_task(self, fw_spec):
        doc = {"name": self["name"],
               "filesystem": env_chk(self.get('filesystem', None), fw_spec),
               "path": self.get("path", os.getcwd())}

        return FWAction(mod_spec=[{'_push': {'vasp_locs': doc}}])


@explicit_serialize
class CopyVaspInputs(FireTaskBase):
    """
    Copy inputs from a previous VASP run directory to the current directory. Additional files can also be specified.

    Note that you must specify either "vasp_dir" or "vasp_loc" of the directory containing the previous VASP run.

    Optional params:
        vasp_dir (str): path to dir (on current filesystem) that contains VASP output files.
        vasp_loc (str OR bool): if True will set most recent vasp_loc. If str search for the most recent vasp_loc with the matching name
        additional_files ([str]): additional files to copy, e.g. ["CHGCAR", "WAVECAR"]. Use $ALL if you just want to copy everything
        contcar_to_poscar(bool): If True, will move CONTCAR to POSCAR

    """

    def run_task(self, fw_spec):
        if 'vasp_dir' not in self and 'vasp_loc' not in self:
            raise ValueError("CopyVaspInputs requires setting either 'vasp_dir' or 'vasp_loc'!")

        vasp_dir = None
        # TODO: separate this function into a utility
        if "vasp_dir" in self:
            vasp_dir = self["vasp_dir"]
        elif self.get("vasp_loc"):
            if isinstance(self["vasp_loc"], basestring):
                for doc in reversed(fw_spec["vasp_locs"]):
                    if doc["name"] == self["vasp_loc_name"]:
                        vasp_dir = doc["path"]
                        break
            else:
                vasp_dir = fw_spec["vasp_locs"][-1]["path"]

        # TODO: add verbosity to logging

        if "$ALL" in self.get("additional_files", []):
            files_to_copy = os.listdir(vasp_dir)
        else:
            files_to_copy = ['INCAR', 'POSCAR', 'KPOINTS', 'POTCAR']
            if self.get("additional_files"):
                files_to_copy.extend(self["additional_files"])

        if self.get("contcar_to_poscar") and "CONTCAR" not in files_to_copy:
            files_to_copy.append("CONTCAR")
            files_to_copy = [f for f in files_to_copy if f != 'POSCAR']  # remove POSCAR

        # TODO: handle gz
        for f in files_to_copy:
            # TODO: handle gz
            dest_fname = 'POSCAR' if f == 'CONTCAR' and self.get("contcar_to_poscar") else f
            prev_path = os.path.join(vasp_dir, f)
            dest_path = os.path.join(os.getcwd(), dest_fname)
            shutil.copy2(prev_path, dest_path)

            # TODO: handle gz
            """
            if '.gz' in dest_file:
                # unzip dest file
                f = gzip.open(dest_file, 'rb')
                file_content = f.read()
                with open(dest_file[0:-3], 'wb') as f_out:
                    f_out.writelines(file_content)
                f.close()
                os.remove(dest_file)
            """