# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os, unittest, glob

from fireworks.core.rocket_launcher import launch_rocket

from atomate.vasp.workflows.base.surface import SurfaceWorkflowCreator
from atomate.vasp.powerups import use_fake_vasp
from atomate.utils.testing import AtomateTest

from pymatgen.core.surface import SlabGenerator, Structure


__author__ = "Richard Tran"
__email__ = "rit001@eng.ucsd.edu"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
ref_dir = os.path.join(module_dir, "..", "..", "test_files")

DEBUG_MODE = False  # If True, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...


class TestSurfaceWorkflow(AtomateTest):

    def setUp(self):
        super(TestSurfaceWorkflow, self).setUp()

        # Set up workflow
        reference_dir = os.path.join(ref_dir, "surface_wf")
        self.ucell = Structure.from_file(os.path.join(reference_dir,
                                                      "Li_mp-135_conventional_unit_cell_k45",
                                                      "inputs", "POSCAR"))
        self.slab = SlabGenerator(self.ucell, (1, 1, 1), 10, 10,
                                  max_normal_search=1).get_slabs()[0]
        self.wfgen = SurfaceWorkflowCreator(db_file=os.path.join(db_dir, "db.json"),
                                            scratch_dir=reference_dir,
                                            run_dir=reference_dir, k_product=45,
                                            vasp_cmd="vasp")
        self.mpid = "mp-135"

    def test_wf(self):

        # Run workflows

        # Test workflow from conventional unit cell
        self.lp.reset(".", require_password=False)
        wf = self.wfgen.from_conventional_unit_cell(self.ucell, 1, naming_tag=self.mpid)
        folder = "Li_mp-135_conventional_unit_cell_k45"
        base_dir = os.path.join(reference_dir, folder)
        wf = use_fake_vasp(wf, {"1": base_dir})
        self.lp.add_wf(wf)
        launch_rocket(self.lp)
        tear_down(base_dir)

        # Test workflow from oriented unit cell
        self.lp.reset(".", require_password=False)
        scale_factor = self.slab.scale_factor
        oriented_ucell = self.slab.oriented_unit_cell
        wf = self.wfgen.from_oriented_unit_cell(oriented_ucell, (1,1,1),
                                                scale_factor, naming_tag=self.mpid)
        folder = "Li_mp-135_bulk_k45_111"
        base_dir = os.path.join(reference_dir, folder)
        wf = use_fake_vasp(wf, {"1": base_dir})
        self.lp.add_wf(wf)
        launch_rocket(self.lp)
        tear_down(base_dir)

        # Test workflow from slab cell
        self.lp.reset(".", require_password=False)
        scale_factor = self.slab.scale_factor
        wf = self.wfgen.from_slab_cell(self.slab, (1,1,1), self.slab.shift,
                                       scale_factor, oriented_ucell, 10,
                                       10, naming_tag=self.mpid)
        folder = "Li_mp-135_slab_k45_s10v10_111_shift0.08333333333333348"
        base_dir = os.path.join(reference_dir, folder)
        wf = use_fake_vasp(wf, {"1": base_dir})
        self.lp.add_wf(wf)
        launch_rocket(self.lp)
        tear_down(base_dir)

        self.lp.reset(".", require_password=False)


def tear_down(base_dir):

    # Remove outputs when done
    os.remove(os.path.join(base_dir, "OUTCAR.relax2.gz"))
    os.remove(os.path.join(base_dir, "CONTCAR.relax2.gz"))
    os.remove(os.path.join(base_dir, "vasprun.xml.relax2.gz"))
    if "FW.json" in glob.glob("*"):
        os.remove(os.path.join(base_dir, "FW.json"))
    if "slab" in base_dir:
        os.remove(os.path.join(base_dir, "LOCPOT.relax2.gz"))



if __name__ == "__main__":
    unittest.main()

