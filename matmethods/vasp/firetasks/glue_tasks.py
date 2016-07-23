# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines tasks that acts as a glue between other vasp firetasks
namely passing the location of current run to the next one and copying files
from previous run directory oto the current one.
"""

import gzip
import os
import re

from pymatgen import MPRester
from pymatgen.io.vasp.sets import get_vasprun_outcar

from fireworks import explicit_serialize, FireTaskBase, FWAction

from matmethods.utils.utils import get_calc_loc, env_chk
from matmethods.utils.fileio import FileClient
from pymatgen.core.structure import Structure
from matmethods.common.firetasks.glue_tasks import GrabFilesFromCalcLoc

__author__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'


@explicit_serialize
class CopyVaspOutputs(FireTaskBase):
    """
    Copy files from a previous VASP run directory to the current directory.
    By default, copies 'INCAR', 'POSCAR', 'KPOINTS', 'POTCAR', 'OUTCAR',
    and 'vasprun.xml'. Additional files, e.g. 'CHGCAR', can also be specified.
    Automatically handles files that have a ".gz" extension (copies and unzips).

    Note that you must specify either "calc_loc" or "calc_dir" to indicate
    the directory containing the previous VASP run.

    Required params:
        (none) - but you must specify either "calc_loc" OR "calc_dir"

    Optional params:
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
            search for the most recent calc_loc with the matching name
        calc_dir (str): path to dir that contains VASP output files.
        filesystem (str): remote filesystem. e.g. username@host
        additional_files ([str]): additional files to copy,
            e.g. ["CHGCAR", "WAVECAR"]. Use $ALL if you just want to copy
            everything
        contcar_to_poscar(bool): If True (default), will move CONTCAR to
            POSCAR (original POSCAR is not copied).
    """

    optional_params = ["calc_loc", "calc_dir", "filesystem", "additional_files", "contcar_to_poscar"]

    def run_task(self, fw_spec):

        if self.get("calc_dir"):  # direct setting of calc dir - no calc_locs or filesystem!
            calc_dir = self["calc_dir"]
            filesystem = None
        elif self.get("calc_loc"):  # search for calc dir and filesystem within calc_locs
            calc_loc = get_calc_loc(self["calc_loc"], fw_spec["calc_locs"])
            calc_dir = calc_loc["path"]
            filesystem = calc_loc["filesystem"]
        else:
            raise ValueError("Must specify either calc_dir or calc_loc!")

        fileclient = FileClient(filesystem=filesystem)
        calc_dir = fileclient.abspath(calc_dir)
        contcar_to_poscar = self.get("contcar_to_poscar", True)

        all_files = fileclient.listdir(calc_dir)

        # determine what files need to be copied
        if "$ALL" in self.get("additional_files", []):
            files_to_copy = all_files
        else:
            files_to_copy = ['INCAR', 'POSCAR', 'KPOINTS', 'POTCAR', 'OUTCAR', 'vasprun.xml']

            if self.get("additional_files"):
                files_to_copy.extend(self["additional_files"])

        if contcar_to_poscar and "CONTCAR" not in files_to_copy:
            files_to_copy.append("CONTCAR")
            files_to_copy = [f for f in files_to_copy if f != 'POSCAR']  # remove POSCAR

        # start file copy
        for f in files_to_copy:
            prev_path_full = os.path.join(calc_dir, f)
            # prev_path = os.path.join(os.path.split(calc_dir)[1], f)
            dest_fname = 'POSCAR' if f == 'CONTCAR' and contcar_to_poscar else f
            dest_path = os.path.join(os.getcwd(), dest_fname)

            relax_ext = ""
            relax_paths = sorted(fileclient.glob(prev_path_full+".relax*"), reverse=True)
            if relax_paths:
                if len(relax_paths) > 9:
                    raise ValueError("CopyVaspOutputs doesn't properly handle >9 relaxations!")
                m = re.search('\.relax\d*', relax_paths[0])
                relax_ext = m.group(0)

            # detect .gz extension if needed - note that monty zpath() did not seem useful here
            gz_ext = ""
            if not (f + relax_ext) in all_files:
                for possible_ext in [".gz", ".GZ"]:
                    if (f + relax_ext + possible_ext) in all_files:
                        gz_ext = possible_ext

            if not (f + relax_ext + gz_ext) in all_files:
                raise ValueError("Cannot find file: {}".format(f))

            # copy the file (minus the relaxation extension)
            fileclient.copy(prev_path_full + relax_ext + gz_ext, dest_path + gz_ext)

            # unzip the .gz if needed
            if gz_ext in ['.gz', ".GZ"]:
                # unzip dest file
                f = gzip.open(dest_path + gz_ext, 'rt')
                file_content = f.read()
                with open(dest_path, 'w') as f_out:
                    f_out.writelines(file_content)
                f.close()
                os.remove(dest_path + gz_ext)

@explicit_serialize
class GetInterpolatedPOSCAR(FireTaskBase):
    """
    Grabs CONTCARS from two previous calculations
    """
    required_params = ["start","end","this_image","nimages"]

    def run_task(self, fw_spec):

        # make folder for poscar interpolation start and end structure files.
        interpolate_folder = '/interpolate'
        if not os.path.exists(os.getcwd()+interpolate_folder):
            os.makedirs(os.getcwd()+interpolate_folder)

            print (os.getcwd()+interpolate_folder)

        # use method of GrabFilesFromCalcLoc to grab files from previous locations.
        GrabFilesFromCalcLoc(calc_dir=None, calc_loc=self.get("start","default"), filenames="CONTCAR",
                             name_prepend="interpolate/", name_append="_0").run_task(fw_spec=fw_spec)
        GrabFilesFromCalcLoc(calc_dir=None, calc_loc=self.get("end","default"), filenames="CONTCAR",
                             name_prepend="interpolate/", name_append="_1").run_task(fw_spec=fw_spec)

        # assuming first calc_dir is polar structure for ferroelectric search

        s1 = Structure.from_file("interpolate/CONTCAR_0")
        s2 = Structure.from_file("interpolate/CONTCAR_1")

        structs = s1.interpolate(s2,self.get('nimages',5),interpolate_lattices=True)

        # save only the interpolation needed for this run
        i = self.get("this_image",0)
        s = structs[i]
        s.to(fmt='POSCAR', filename=os.getcwd()+"/POSCAR")

        # will want to call this method after getting first VASP set

@explicit_serialize
class CheckStability(FireTaskBase):
    """
    Checks the stability of the entry against the Materials Project database.
    If the stability is less than the cutoff (default is 0.1 eV/atom), then
    the task will return a FWAction that will defuse all remaining tasks.

    Required params:
        (none) - but your MAPI key must be set as an environ var in this case

    Optional params:
        ehull_cutoff: (float) energy in eV/atom to use as ehull cutoff. Default
            is 0.05 eV/atom.
        MAPI_KEY: (str) set MAPI key directly. Supports env_chk.
        calc_dir: (str) string to path containing vasprun.xml (default currdir)
    """

    required_params = []
    optional_params = ["ehull_cutoff", "MAPI_KEY", "calc_dir"]

    def run_task(self, fw_spec):

        mpr = MPRester(env_chk(self.get("MAPI_KEY"), fw_spec))
        vasprun, outcar = get_vasprun_outcar(self.get("calc_dir", "."), parse_dos=False, parse_eigen=False)

        my_entry = vasprun.get_computed_entry(inc_structure=False)
        stored_data = mpr.get_stability([my_entry])[0]

        if stored_data["e_above_hull"] > self.get("ehull_cutoff", 0.05):
            return FWAction(stored_data=stored_data, exit=True, defuse_workflow=True)

        else:
            return FWAction(stored_data=stored_data)
