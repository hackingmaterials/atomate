# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
This module defines tasks that acts as a glue between other vasp firetasks
namely passing the location of current run to the next one and copying files
from previous run directory oto the current one.
"""

import gzip
import os
import shutil
import re
from glob import glob

from fireworks import explicit_serialize, FireTaskBase, FWAction

from matmethods.utils.utils import env_chk
from matmethods.vasp.vasp_utils import get_vasp_dir

__author__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'


@explicit_serialize
class PassVaspLocs(FireTaskBase):
    """
    Passes the vasp_locs key. Should be called in the same FireWork as a
    VASP run. This passes information about where the current run is located
    for the next FireWork.

    Required params:
        name: descriptive name for this VASP file/dir

    Optional params:
        filesystem: name of filesystem. Supports env_chk. defaults to None
        path: The path to the directory containing the VASP run. defaults to
            current working directory.
    """

    required_params = ["name"]
    optional_params = ["filesystem", "path"]

    def run_task(self, fw_spec):
        doc = {"name": self["name"],
               "filesystem": env_chk(self.get('filesystem', None), fw_spec),
               "path": self.get("path", os.getcwd())}

        return FWAction(mod_spec=[{'_push': {'vasp_locs': doc}}])


@explicit_serialize
class CopyVaspOutputs(FireTaskBase):
    """
    Copy outputs from a previous VASP run directory to the current directory.
    Additional files, e.g. CHGCAR, can also be specified.
    Automatically handles files that have a ".gz" extension (copies and unzips).

    Note that you must specify either "vasp_loc" or "vasp_dir" to indicate
    the directory containing the previous VASP run.

    Required params:
        (none) - but you must specify either "vasp_loc" OR "vasp_dir"

    Optional params:
        vasp_loc (str OR bool): if True will set most recent vasp_loc. If str
            search for the most recent vasp_loc with the matching name
        vasp_dir (str): path to dir (on current filesystem) that contains VASP
            output files.
        additional_files ([str]): additional files to copy,
            e.g. ["CHGCAR", "WAVECAR"]. Use $ALL if you just want to copy
            everything
        contcar_to_poscar(bool): If True (default), will move CONTCAR to
            POSCAR (original POSCAR is not copied).
    """

    optional_params = ["vasp_loc", "vasp_dir", "additional_files",
                       "contcar_to_poscar"]

    def run_task(self, fw_spec):

        vasp_dir = get_vasp_dir(self, fw_spec)
        contcar_to_poscar = self.get("contcar_to_poscar", True)

        # determine what files need to be copied
        if "$ALL" in self.get("additional_files", []):
            files_to_copy = os.listdir(vasp_dir)
        else:
            files_to_copy = ['INCAR', 'POSCAR', 'KPOINTS', 'POTCAR', 'OUTCAR',
                             'vasprun.xml']
            if self.get("additional_files"):
                files_to_copy.extend(self["additional_files"])

        if contcar_to_poscar and "CONTCAR" not in files_to_copy:
            files_to_copy.append("CONTCAR")
            files_to_copy = [f for f in files_to_copy if
                             f != 'POSCAR']  # remove POSCAR

        # start file copy
        for f in files_to_copy:
            prev_path = os.path.join(vasp_dir, f)
            dest_fname = 'POSCAR' if f == 'CONTCAR' and contcar_to_poscar else f
            dest_path = os.path.join(os.getcwd(), dest_fname)

            # detect .relax## if needed - uses last relaxation (up to 9 relaxations)
            relax_ext = ""
            # TODO: modify sorting for >9 relaxations
            relax_path = sorted(glob(prev_path+".relax*"), reverse=True)
            if relax_path:
                # TODO: change \d to \d* for >9 relaxations
                m = re.search('\.relax\d', relax_path[0])
                relax_ext = m.group(0)

            # detect .gz extension if needed - note that monty zpath() did not
            # seem useful here
            gz_ext = ""
            if not os.path.exists(prev_path + relax_ext):
                for possible_ext in [".gz", ".GZ"]:
                    if os.path.exists(prev_path + relax_ext + possible_ext):
                        gz_ext = possible_ext

            if not os.path.exists(prev_path + relax_ext + gz_ext):
                raise ValueError("Cannot find file: {}".format(prev_path))

            # copy the file (minus the relaxation extension)
            shutil.copy2(prev_path + relax_ext + gz_ext, dest_path + gz_ext)

            # unzip the .gz if needed
            if gz_ext == '.gz' or gz_ext == ".GZ":
                # unzip dest file
                f = gzip.open(dest_path + gz_ext, 'rb')
                file_content = f.read()
                with open(dest_path, 'wb') as f_out:
                    f_out.writelines(file_content)
                f.close()
                os.remove(dest_path + gz_ext)
