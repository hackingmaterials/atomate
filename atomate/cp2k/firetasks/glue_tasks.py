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

from pymatgen.io.cp2k.inputs import Cp2kInput, Coord, Cell
from pymatgen.io.cp2k.outputs import parse_structures

from fireworks import explicit_serialize, FiretaskBase, FWAction

from atomate.utils.utils import env_chk, get_logger
from atomate.common.firetasks.glue_tasks import get_calc_loc, PassResult, \
    CopyFiles, CopyFilesFromCalcLoc

logger = get_logger(__name__)

__author__ = 'Nicholas Winner'
__email__ = 'nwinner@berkeley.edu'


@explicit_serialize
class PassStructureToCp2kSet:

    required_params = ['trajectory_file', 'lattice_file']

    def run_task(self, fw_spec):
        cp2k_input_set = fw_spec.get('cp2k_input_set')
        ci = Cp2kInput.from_dict(cp2k_input_set)

        trajectory_file = self.get('trajectory_file')
        lattice_file = self.get('lattice_file')
        final = parse_structures(trajectory_file=trajectory_file, lattice_file=lattice_file, final=True)
        ci['FORCE_EVAL']['SUBSYS']['COORD'] = Coord(final)
        ci['FORCE_EVAL']['SUBSYS']['CELL'] = Cell(final.lattice)

        fw_spec['cp2k_input_set'] = ci.as_dict()
