# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
import shutil

from atomate.qchem.firetasks.run_calc import RunQChemCustodian
from atomate.qchem.firetasks.write_inputs import WriteInputFromIOSet
from atomate.qchem.firetasks.parse_outputs import QChemToDb
from atomate.qchem.fireworks.core import FrequencyFlatteningOptimizeFW
from atomate.qchem.workflows.base.double_FF_opt import get_wf_double_FF_opt
from fireworks import Firework, Workflow, FWorker
from fireworks.core.rocket_launcher import rapidfire
from atomate.utils.testing import AtomateTest
from pymatgen.core import Molecule
from pymatgen.io.qchem_io.outputs import QCOutput
from pymatgen.io.qchem_io.inputs import QCInput
import numpy as np

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")

class TestDoubleFFOpt(AtomateTest):
    @classmethod
    def setUpClass(cls):
        out_file = os.path.join(module_dir, "..", "..", "test_files", "FF_working", "test.qout.opt_0")
        qc_out = QCOutput(filename=out_file)
        cls.act_mol = qc_out.data["initial_molecule"]

    def tearDown(self):
    	pass

    # def tearDown(self):
    #     # this removes the scratch dir made by AtomateTest
    #     shutil.rmtree(self.scratch_dir)
    #     # this removes the file that gets written
    #     for x in ["mol.qin"]:
    #         if os.path.exists(os.path.join(module_dir, x)):
    #             os.remove(os.path.join(module_dir, x))

    def test_double_FF_opt(self):
    	wf = get_wf_double_FF_opt(mymol, 10.0)
    	self.lp.add_wf(wf)
    	rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))