# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

from fireworks import FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.vasp.powerups import use_fake_vasp
from atomate.vasp.workflows.base.surface import get_wflow_from_mpid, get_fw_from_ucell
from atomate.utils.testing import AtomateTest

from pymatgen import Structure, Molecule, Lattice
from pymatgen.core.surface import generate_all_slabs

__author__ = "Richard Tran"
__email__ = "rit001@eng.ucsd.edu"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
reference_dir = os.path.join(module_dir, "..", "..", "test_files", "surface_wf")

DEBUG_MODE = False  # If True, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...
_write_task_docs = False # Test developer option: defaults to False, need to be True only once


class TestSurfaceWorkflow(AtomateTest):

    def setUp(self):
        super(TestSurfaceWorkflow, self).setUp()

        self.struct_ir = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3.875728),
                                                   ["Ir"], [[0, 0, 0]])
        self.wf = get_fw_from_ucell(self.struct_ir, "vasp", ".",
                                    os.path.join(db_dir, db.json),
                                    inc_conv_ucell=True)
        self.name = "%s_%s_conventional_k%s" %("Ni", "--", 50)
        self.cwd = os.getcwd()

    def run_wf(self):
        ref_dir = os.path.join(self.cwd, self.name)
        use_fake_vasp(original_wf, {"Ni_--": ref_dir})



if __name__ == "__main__":
    unittest.main()
