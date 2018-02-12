# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import
#
import os
import unittest
#
# from fireworks import FWorker
# from fireworks.core.rocket_launcher import rapidfire
#
# from atomate.vasp.powerups import use_fake_vasp
# # from atomate.vasp.workflows.base.surface import get_wflow_from_mpid, get_fw_from_ucell
from atomate.utils.testing import AtomateTest
# from atomate.vasp.fireworks.core import SurfCalcFW
#
from pymatgen import Structure, Molecule, Lattice
# from pymatgen.core.surface import generate_all_slabs

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
        self.struct_li = Structure.from_spacegroup("Im-3m", Lattice.cubic(3.478000),
                                                   ["Li"], [[0, 0, 0]])
        # self.wf = get_fw_from_ucell(self.struct_ir, "vasp", ".",
        #                             os.path.join(db_dir, db.json),
        #                             inc_conv_ucell=True)
        # self.name = "%s_%s_conventional_k%s" %("Ni", "--", 50)
        # self.cwd = os.getcwd()

    def test_wf(self):
        # ref_dir = os.path.join(self.cwd, self.name)
        fw = SurfCalcFW(self.struct_li, "conventional_unit_cell", "mp-135",
                        k_product=50, miller_index=None, min_slab_size=None,
                        min_vac_size=None, reconstruction_name=None, shift=None,
                        user_incar_settings={}, ediffg=-0.02, vasp_cmd="vasp")
        # wf = Workflow([fw], name="test")
        #
        # wf = use_fake_vasp(wf, {"Li-mp-135_conventional_unit_cell_k50": ref_dir})
        # wf = add_namefile(wf)
        # self.lp.add_wf(wf)
        # rapidfire(self.lp)


if __name__ == "__main__":
    unittest.main()

    # def cleanup(self):
    #
    #     pass
    #
    # def simulate_vasprun(self):
    #
    #     pass



