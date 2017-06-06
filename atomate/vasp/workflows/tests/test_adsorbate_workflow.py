# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import shutil
import unittest

from fireworks import LaunchPad, FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.vasp.powerups import use_fake_vasp
from atomate.vasp.workflows.base.adsorption import get_wf_surface
from atomate.utils.testing import AtomateTest

from pymatgen import SETTINGS, Structure, Molecule, Lattice
from pymatgen.core.surface import generate_all_slabs

__author__ = 'Kiran Mathew, Joseph Montoya'
__email__ = 'montoyjh@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
ref_dir = os.path.join(module_dir, "..", "..", "test_files")

DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...


class TestAdsorptionWorkflow(AtomateTest):
    @classmethod
    def setUpClass(cls):
        if not SETTINGS.get("PMG_VASP_PSP_DIR"):
            SETTINGS["PMG_VASP_PSP_DIR"] = os.path.join(module_dir, "..", "..", "tests", "..", "..", "test_files")
            print('This system is not set up to run VASP jobs. '
                  'Please set PMG_VASP_PSP_DIR variable in your ~/.pmgrc.yaml file.')

        cls.struct_ir = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3.875728), ["Ir"], [[0, 0, 0]])
        cls.scratch_dir = os.path.join(module_dir, "scratch")
        sgp = {"max_index": 1, "min_slab_size": 7.0, "min_vacuum_size": 20.0}
        slabs = generate_all_slabs(cls.struct_ir, **sgp)
        slabs = [slab for slab in slabs if slab.miller_index==(1, 0, 0)]
        sgp.pop("max_index")
        cls.wf_1 = get_wf_surface(slabs, [Molecule("H", [[0, 0, 0]])], cls.struct_ir, sgp,
                                  db_file=os.path.join(db_dir, "db.json"))

    def setUp(self):
        if os.path.exists(self.scratch_dir):
            shutil.rmtree(self.scratch_dir)
        os.makedirs(self.scratch_dir)
        os.chdir(self.scratch_dir)
        try:
            self.lp = LaunchPad.from_file(os.path.join(db_dir, "my_launchpad.yaml"))
            self.lp.reset("", require_password=False)
        except:
            raise unittest.SkipTest(
                'Cannot connect to MongoDB! Is the database server running? '
                'Are the credentials correct?')

    def tearDown(self):
        if not DEBUG_MODE:
            shutil.rmtree(self.scratch_dir)
            self.lp.reset("", require_password=False)
            db = self.get_task_database()
            for coll in db.collection_names():
                if coll != "system.indexes":
                    db[coll].drop()
            os.chdir(module_dir)

    def _simulate_vasprun(self, wf):
        reference_dir = os.path.abspath(os.path.join(ref_dir, "adsorbate_wf"))
        ir_ref_dirs = {"Ir-structure optimization": os.path.join(reference_dir, "1"),
                       "Ir-Ir_(1, 0, 0) slab optimization": os.path.join(reference_dir, "2"),
                       "Ir-H1-Ir_(1, 0, 0) adsorbate optimization 0": os.path.join(reference_dir, "3"),
                       "Ir-H1-Ir_(1, 0, 0) adsorbate optimization 1": os.path.join(reference_dir, "4"),
                       "Ir-H1-Ir_(1, 0, 0) adsorbate optimization 2": os.path.join(reference_dir, "5")}
        return use_fake_vasp(wf, ir_ref_dirs, params_to_check=["ENCUT", "ISIF", "IBRION"])

    def _check_run(self, d, mode):
        if mode not in ["H1-Ir_(1, 0, 0) adsorbate optimization 1"]:
            raise ValueError("Invalid mode!")

        if "adsorbate" in mode:
            self.assertEqual(d["formula_reduced_abc"], "H1 Ir16")
        # Check relaxation of adsorbate
        # Check slab calculations
        # Check structure optimization

    def test_wf(self):
        wf = self._simulate_vasprun(self.wf_1)

        self.assertEqual(len(self.wf_1.fws), 5)
        # check vasp parameters for ionic relaxation
        ads_vis = [fw.tasks[1]['vasp_input_set']
                   for fw in self.wf_1.fws if "adsorbate" in fw.name]
        assert all([vis.incar['EDIFFG']==-0.01 for vis in ads_vis])
        assert all([vis.incar['ISIF']==2 for vis in ads_vis])
        self.lp.add_wf(wf)
        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))

        # check relaxation
        d = self.get_task_collection().find_one({"task_label": "H1-Ir_(1, 0, 0) adsorbate optimization 1"})
        self._check_run(d, mode="H1-Ir_(1, 0, 0) adsorbate optimization 1")

        wf = self.lp.get_wf_by_fw_id(1)
        self.assertTrue(all([s == 'COMPLETED' for s in wf.fw_states.values()]))


if __name__ == "__main__":
    unittest.main()
