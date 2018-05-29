# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
import shutil

from atomate.qchem.firetasks.geo_transformations import RotateTorsion
from atomate.qchem.firetasks.write_inputs import WriteInputFromIOSet
from atomate.qchem.firetasks.parse_outputs import QChemToDb
from fireworks import Firework, Workflow, FWorker
from fireworks.core.rocket_launcher import rapidfire
from atomate.utils.testing import AtomateTest
from pymatgen.core import Molecule
from pymatgen.io.qchem_io.outputs import QCOutput
from pymatgen.io.qchem_io.inputs import QCInput
import numpy as np

__author__ = 'Brandon Wood'
__email__ = 'b.wood@berkeley.edu'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")

DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test


class TestTorsionPotential(AtomateTest):
    @classmethod
    def setUpClass(cls):
        # define starting molecule and torsion potential workflow object

        out_file = os.path.join(module_dir, "..", "..", "test_files", "FF_working", "test.qout.opt_1")
        qc_out = QCOutput(filename=out_file)
        cls.act_mol = qc_out.data["molecule_from_optimized_geometry"]
