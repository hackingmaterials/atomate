# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import json
import os
import shutil
import unittest

import numpy as np

from pymongo import MongoClient

from monty.os.path import which

from atomate.vasp.workflows.base.magnetism import MagneticOrderingsWF
from atomate.utils.testing import AtomateTest
from atomate.vasp.powerups import use_fake_vasp

from pymatgen import Structure

__author__ = 'Matthew Horton'
__email__ = 'mkhorton@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
ref_dir = os.path.join(module_dir, "..", "..", "test_files", "magnetism_wf")

enum_cmd = which('enum.x') or which('multienum.x')
makestr_cmd = which('makestr.x') or which('makeStr.x') or which('makeStr.py')
enumlib_present = enum_cmd and makestr_cmd


@unittest.skipIf(not enumlib_present, "enumlib not present")
class TestMagnetismWorkflow(unittest.TestCase):

    def test_magnetic_orderings_workflow(self):

        # simple afm
        structure = Structure.from_file(os.path.join(ref_dir, "LaMnO3.json"))
        wf = MagneticOrderingsWF(structure)
        self.assertEqual(wf.input_origin, "afm")

        # ferrimagnetic (Cr produces net spin)
        structure = Structure.from_file(os.path.join(ref_dir, "Cr2NiO4.json"))
        wf = MagneticOrderingsWF(structure)
        self.assertEqual(wf.input_origin, "ferri_by_Cr")

        # antiferromagnetic on single magnetic site
        structure = Structure.from_file(os.path.join(ref_dir, "Cr2WO6.json"))
        wf = MagneticOrderingsWF(structure)
        self.assertEqual(wf.input_origin, "afm_by_Cr")

        # afm requiring large cell size (enable for development, too slow for CI)
        #structure = Structure.from_file(os.path.join(ref_dir, "CuO.json"))
        #wf = MagneticOrderingsWF(structure, default_magmoms={'Cu': 1.73},
        #                         transformation_kwargs={'max_cell_size': 4})
        #self.assertEqual(wf.input_origin, "afm")

        # antiferromagnetic by structural motif
        structure = Structure.from_file(os.path.join(ref_dir, "Ca3Co2O6.json"))
        wf = MagneticOrderingsWF(structure, strategies=('antiferromagnetic_by_motif', ),
                                 # this example just misses default cut-off, so do not truncate
                                 truncate_by_symmetry=False,
                                 transformation_kwargs={'max_cell_size': 2})
        self.assertEqual(wf.input_origin, "afm_by_motif_2a")
