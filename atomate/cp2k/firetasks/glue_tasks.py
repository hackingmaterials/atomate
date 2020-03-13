# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import glob

from pymatgen.analysis.elasticity.strain import Strain
from pymatgen.io.vasp import Vasprun, zpath

"""
This module defines tasks that acts as a glue between other vasp Firetasks to allow communication
between different Firetasks and Fireworks. This module also contains tasks that affect the control
flow of the workflow, e.g. tasks to check stability or the gap is within a certain range.
"""

import gzip
import os
import re
import pickle
import glob
from monty.io import zopen

from pymatgen.io.cp2k.sets import Cp2kInputSet
from pymatgen.io.cp2k.inputs import Cp2kInput, Coord, Cell
from pymatgen.io.cp2k.outputs import Cp2kOutput

from fireworks import explicit_serialize, FiretaskBase, FWAction

from atomate.utils.utils import env_chk, get_logger
from atomate.common.firetasks.glue_tasks import get_calc_loc, PassResult, \
    CopyFiles, CopyFilesFromCalcLoc

logger = get_logger(__name__)

__author__ = 'Nicholas Winner'
__email__ = 'nwinner@berkeley.edu'


@explicit_serialize
class UpdateStructureFromPrevCalc(FiretaskBase):
    """
    Using the location of a previous calculation. The CP2K output parser will
    get the final structure from a previous calculation and update this FW's
    cp2k_input_set using it.
    """

    required_params = ['prev_calc_loc']
    optional_params = ['cp2k_output_file']

    def run_task(self, fw_spec):
        calc_loc = get_calc_loc(self.get('prev_calc_loc'),
                                fw_spec["calc_locs"]) if self.get(
            "prev_calc_loc") else {}

        cp2k_input_set = fw_spec.get('cp2k_input_set')
        ci = Cp2kInputSet.from_dict(cp2k_input_set)

        if self.get('cp2k_output_file'):
            out = Cp2kOutput(os.path.join(calc_loc['path'], self.get('cp2k_output_file')))
        else:
            out = Cp2kOutput(glob.glob(calc_loc['path']+'/cp2k.out*')[0])

        out.parse_structures()
        ci['FORCE_EVAL']['SUBSYS']['COORD'] = Coord(out.final_structure)
        ci['FORCE_EVAL']['SUBSYS']['CELL'] = Cell(out.final_structure.lattice)
        fw_spec['cp2k_input_set'] = ci.as_dict()


@explicit_serialize
class CopyCp2kOutputs(CopyFiles):
    """
    Copy CP2K outputs from from one location to another, unzipping them if necessary.
    Unlike VASP, which might have CONTCAR copied to POSCAR in order to continue a calculation,
    Cp2k doesn't exactly use this file system. What will generally be used for is to copy
    a .wfn file in order to restart a calculation or as an initial guess for a hybrid calculation.


    Note that you must specify either "calc_loc" or "calc_dir" to indicate
    the directory containing the previous VASP run.

    Required params:
        (none) - but you must specify either "calc_loc" OR "calc_dir"

    Optional params:
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
            search for the most recent calc_loc with the matching name
        calc_dir (str): path to dir that contains VASP output files.
        filesystem (str): remote filesystem. e.g. username@host
    """

    optional_params = ["files_to_copy", "calc_loc", "calc_dir", "filesystem"]

    def run_task(self, fw_spec):

        calc_loc = get_calc_loc(self["calc_loc"],
                                fw_spec["calc_locs"]) if self.get(
            "calc_loc") else {}

        files_to_copy = self.get('files_to_copy', [])
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
            dest_fname = f
            dest_path = os.path.join(self.to_dir, dest_fname)

            # detect .gz extension if needed - note that monty zpath() did not seem useful here
            gz_ext = ""
            if not f in all_files:
                for possible_ext in [".gz", ".GZ"]:
                    if (f + possible_ext) in all_files:
                        gz_ext = possible_ext

            if not (f + gz_ext) in all_files:
                raise ValueError("Cannot find file: {}".format(f))

            # copy the file (minus the relaxation extension)
            self.fileclient.copy(prev_path_full + gz_ext,
                                 dest_path + gz_ext)

            # unzip the .gz if needed
            if gz_ext in ['.gz', ".GZ"]:
                # unzip dest file
                try:
                    f = zopen(dest_path + gz_ext, 'rt')
                    file_content = f.read()
                except (UnicodeDecodeError, AttributeError):
                    f = zopen(dest_path + gz_ext, 'rb')
                    file_content = f.read()
                if isinstance(file_content, (bytes, bytearray)):
                    with open(dest_path, 'wb') as f_out:
                        f_out.write(file_content)
                else:
                    with open(dest_path, 'w') as f_out:
                        f_out.write(file_content)

                f.close()
                os.remove(dest_path + gz_ext)
