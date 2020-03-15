# coding: utf-8

import glob

from pymatgen.analysis.elasticity.strain import Strain
from pymatgen.io.vasp import Vasprun
from monty.os.path import zpath

"""
This module defines tasks that acts as a glue between other vasp Firetasks to allow communication
between different Firetasks and Fireworks. This module also contains tasks that affect the control
flow of the workflow, e.g. tasks to check stability or the gap is within a certain range.
"""

import shutil
import gzip
import os
import re

from pymatgen import MPRester
from pymatgen.io.vasp.sets import get_vasprun_outcar
from pymatgen.core.structure import Structure

from fireworks import explicit_serialize, FiretaskBase, FWAction

from atomate.utils.utils import env_chk, get_logger
from atomate.common.firetasks.glue_tasks import get_calc_loc, PassResult, \
    CopyFiles, CopyFilesFromCalcLoc

logger = get_logger(__name__)

__author__ = 'Anubhav Jain, Kiran Mathew'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov'


@explicit_serialize
class CopyVaspOutputs(CopyFiles):
    """
    Copy files from a previous VASP run directory to the current directory.
    By default, copies 'INCAR', 'POSCAR' (default: via 'CONTCAR'), 'KPOINTS',
    'POTCAR', 'OUTCAR', and 'vasprun.xml'. Additional files, e.g. 'CHGCAR',
    can also be specified. Automatically handles files that have a ".gz"
    extension (copies and unzips).

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

    optional_params = ["calc_loc", "calc_dir", "filesystem", "additional_files",
                       "contcar_to_poscar"]

    def run_task(self, fw_spec):

        calc_loc = get_calc_loc(self["calc_loc"],
                                fw_spec["calc_locs"]) if self.get(
            "calc_loc") else {}

        # determine what files need to be copied
        files_to_copy = None
        if not "$ALL" in self.get("additional_files", []):
            files_to_copy = ['INCAR', 'POSCAR', 'KPOINTS', 'POTCAR', 'OUTCAR',
                             'vasprun.xml']
            if self.get("additional_files"):
                files_to_copy.extend(self["additional_files"])

        # decide between poscar and contcar
        contcar_to_poscar = self.get("contcar_to_poscar", True)
        if contcar_to_poscar and "CONTCAR" not in files_to_copy:
            files_to_copy.append("CONTCAR")
            files_to_copy = [f for f in files_to_copy if
                             f != 'POSCAR']  # remove POSCAR

        # setup the copy
        self.setup_copy(self.get("calc_dir", None),
                        filesystem=self.get("filesystem", None),
                        files_to_copy=files_to_copy, from_path_dict=calc_loc)
        # do the copying
        self.copy_files()

    def copy_files(self):
        all_files = self.fileclient.listdir(self.from_dir)
        # start file copy
        for f in self.files_to_copy:
            prev_path_full = os.path.join(self.from_dir, f)
            dest_fname = 'POSCAR' if f == 'CONTCAR' and self.get(
                "contcar_to_poscar", True) else f
            dest_path = os.path.join(self.to_dir, dest_fname)

            relax_ext = ""
            relax_paths = sorted(
                self.fileclient.glob(prev_path_full + ".relax*"))
            if relax_paths:
                if len(relax_paths) > 9:
                    raise ValueError(
                        "CopyVaspOutputs doesn't properly handle >9 relaxations!")
                m = re.search('\.relax\d*', relax_paths[-1])
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
            self.fileclient.copy(prev_path_full + relax_ext + gz_ext,
                                 dest_path + gz_ext)

            # unzip the .gz if needed
            if gz_ext in ['.gz', ".GZ"]:
                # unzip dest file
                with open(dest_path, 'wb') as f_out, gzip.open(dest_path + gz_ext, 'rb') as f_in:
                    shutil.copyfileobj(f_in, f_out)
                os.remove(dest_path + gz_ext)


@explicit_serialize
class CheckStability(FiretaskBase):
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
        vasprun, outcar = get_vasprun_outcar(self.get("calc_dir", "."),
                                             parse_dos=False,
                                             parse_eigen=False)

        my_entry = vasprun.get_computed_entry(inc_structure=False)
        stored_data = mpr.get_stability([my_entry])[0]

        if stored_data["e_above_hull"] > self.get("ehull_cutoff", 0.05):
            logger.info("CheckStability: failed test!")
            return FWAction(stored_data=stored_data, exit=True,
                            defuse_workflow=True)

        else:
            return FWAction(stored_data=stored_data)


@explicit_serialize
class CheckBandgap(FiretaskBase):
    """
    Checks the band gap of an entry. If band gap is >min_gap or <max_gap, then
    the task will return a FWAction that will defuse all remaining tasks.

    Required params:
        (none) - but you should set either min_gap or max_gap

    Optional params:
        min_gap: (float) minimum gap energy in eV to proceed
        max_gap: (float) maximum gap energy in eV to proceed
        vasprun_path: (str) path to vasprun.xml file
    """

    required_params = []
    optional_params = ["min_gap", "max_gap", "vasprun_path"]

    def run_task(self, fw_spec):
        vr_path = zpath(self.get("vasprun_path", "vasprun.xml"))
        min_gap = self.get("min_gap", None)
        max_gap = self.get("max_gap", None)

        if not os.path.exists(vr_path):
            relax_paths = sorted(glob.glob(vr_path + ".relax*"))
            if relax_paths:
                if len(relax_paths) > 9:
                    raise ValueError(
                        "CheckBandgap doesn't properly handle >9 relaxations!")
                vr_path = relax_paths[-1]

        logger.info("Checking the gap of file: {}".format(vr_path))
        vr = Vasprun(vr_path, parse_potcar_file=False)
        gap = vr.get_band_structure().get_band_gap()["energy"]
        stored_data = {"band_gap": gap}
        logger.info(
            "The gap is: {}. Min gap: {}. Max gap: {}".format(gap, min_gap,
                                                              max_gap))

        if (min_gap and gap < min_gap) or (max_gap and gap > max_gap):
            logger.info("CheckBandgap: failed test!")
            return FWAction(stored_data=stored_data, exit=True,
                            defuse_workflow=True)

        return FWAction(stored_data=stored_data)


@explicit_serialize
class GetInterpolatedPOSCAR(FiretaskBase):
    """
    Grabs CONTCARS from two previous calculations to create interpolated
    structure.

    The code gets the CONTCAR locations using get_calc_loc of two calculations
    indicated by the start and end params, creates a folder named "interpolate"
    in the current FireWork directory, and copies the two CONTCARs to this folder.
    The two CONTCARs are then used to create nimages interpolated structures using
    pymatgen.core.structure.Structure.interpolate. Finally, the structure indicated
    by this_image is written as a POSCAR file.

    Required params:
        start (str): name of fw for start of interpolation.
        end (str): name of fw for end of interpolation.
        this_image (int): which interpolation this is.
        nimages (int) : number of interpolations.

    Optional params:
        autosort_tol (float): parameter used by Structure.interpolate.
          a distance tolerance in angstrom in which to automatically
          sort end_structure to match to the closest
          points in this particular structure. Default is 0.0.

    """
    required_params = ["start", "end", "this_image", "nimages"]
    optional_params = ["autosort_tol"]

    def run_task(self, fw_spec):
        structure = self.interpolate_poscar(fw_spec)
        structure.to(fmt="POSCAR", filename=os.path.join(os.getcwd(), "POSCAR"))

    def interpolate_poscar(self, fw_spec):
        # make folder for poscar interpolation start and end structure files.
        interpolate_folder = 'interpolate'
        if not os.path.exists(os.path.join(os.getcwd(), interpolate_folder)):
            os.makedirs(os.path.join(os.getcwd(), interpolate_folder))

        # use CopyFilesFromCalcLoc to get files from previous locations.
        CopyFilesFromCalcLoc(calc_loc=self["start"],
                             filenames=["CONTCAR"],
                             name_prepend=interpolate_folder + os.sep,
                             name_append="_0").run_task(fw_spec=fw_spec)
        CopyFilesFromCalcLoc(calc_loc=self["end"],
                             filenames=["CONTCAR"],
                             name_prepend=interpolate_folder + os.sep,
                             name_append="_1").run_task(fw_spec=fw_spec)

        # assuming first calc_dir is polar structure for ferroelectric search
        s1 = Structure.from_file(os.path.join(interpolate_folder, "CONTCAR_0"))
        s2 = Structure.from_file(os.path.join(interpolate_folder, "CONTCAR_1"))

        structs = s1.interpolate(s2, self["nimages"], interpolate_lattices=True,
                                 autosort_tol=self.get("autosort_tol", 0.0))

        # save only the interpolation needed for this run

        i = self.get("this_image")
        return structs[i]


def pass_vasp_result(pass_dict=None, calc_dir='.', filename="vasprun.xml.gz",
                     parse_eigen=False,
                     parse_dos=False, **kwargs):
    """
    Function that gets a PassResult firework corresponding to output from a Vasprun.  Covers
    most use cases in which user needs to pass results from a vasp run to child FWs
    (e. g. analysis FWs)

    pass_vasp_result(pass_dict={'stress': ">>ionic_steps.-1.stress"})

    Args:
        pass_dict (dict): dictionary designating keys and values to pass
            to child fireworks.  If value is a string beginning with '>>',
            the firework will search the parsed VASP output dictionary
            for the designated property by following the sequence of keys
            separated with periods, e. g. ">>ionic_steps.-1.stress" is used
            to designate the stress from the last ionic_step. If the value
            is not a string or does not begin with ">>" or "a>>" (for an
            object attribute, rather than nested key of .as_dict() conversion),
            it is passed as is.  Defaults to pass the computed entry of
            the Vasprun.
        calc_dir (str): path to dir that contains VASP output files, defaults
            to '.', e. g. current directory
        filename (str): filename for vasp xml file to parse, defaults to
            "vasprun.xml.gz"
        parse_eigen (bool): flag on whether or not to parse eigenvalues,
            defaults to false
        parse_eigen (bool): flag on whether or not to parse dos,
            defaults to false
        **kwargs (keyword args): other keyword arguments passed to PassResult
            e.g. mod_spec_key or mod_spec_cmd

    """
    pass_dict = pass_dict or {"computed_entry": "a>>get_computed_entry"}
    parse_kwargs = {"filename": filename, "parse_eigen": parse_eigen,
                    "parse_dos": parse_dos}
    return PassResult(pass_dict=pass_dict, calc_dir=calc_dir,
                      parse_kwargs=parse_kwargs,
                      parse_class="pymatgen.io.vasp.outputs.Vasprun", **kwargs)
