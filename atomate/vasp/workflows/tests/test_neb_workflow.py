# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import yaml
import os
import shutil
import unittest
import copy

from fireworks.core.launchpad import LaunchPad
from fireworks.core.fworker import FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.vasp.powerups import use_fake_vasp
from atomate.vasp.workflows.presets.core import wf_nudged_elastic_band

from pymatgen import SETTINGS
from pymatgen.core import Structure
from pymatgen.util.testing import PymatgenTest

try:
    from pymatgen_diffusion.neb.io import get_endpoints_from_index
    pmgd = True
except ImportError:
    pmgd = False

__author__ = "Hanmei Tang, Iek-Heng Chu"
__email__ = 'hat003@eng.ucsd.edu, ihchu@eng.ucsd.edu'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
ref_dir = os.path.join(module_dir, "..", "..", "test_files")

DEBUG_MODE = False
LAUNCHPAD_RESET = True


@unittest.skipIf(not pmgd, "pymatgen-diffusion not installed, so skipping...")
class TestNudgedElasticBandWorkflow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        1) Basic check for pymatgen configurations.
        2) Setup all test workflow.
        """
        if not SETTINGS.get("PMG_VASP_PSP_DIR"):
            SETTINGS["PMG_VASP_PSP_DIR"] = os.path.join(module_dir, "..", "..", "tests",
                                                        "..", "..", "test_files")
            print('This system is not set up to run VASP jobs. '
                  'Please set PMG_VASP_PSP_DIR variable in '
                  'your ~/.pmgrc.yaml file.')

        # Structures used for test:
        parent = PymatgenTest.get_structure("Li2O")
        parent.remove_oxidation_states()
        parent.make_supercell(2)
        ep0, ep1 = get_endpoints_from_index(parent, [0, 1])
        neb_dir = [os.path.join(module_dir, "..", "..", "test_files", "neb_wf", "4", "inputs", "{:02}",
                                "POSCAR").format(i) for i in range(5)]
        cls.structures = [Structure.from_file(n) for n in neb_dir]
        cls.scratch_dir = os.path.join(module_dir, "scratch")

        # Run fake vasp
        test_yaml = "../../test_files/neb_wf/config/neb_unittest.yaml"
        with open(test_yaml, 'r') as stream:
            cls.config = yaml.load(stream)
            # Use scratch directory as destination directory for testing
            cls.config["common_params"]["_fw_env"] = {"run_dest_root": cls.scratch_dir}

        # Config 1: The parent structure & two endpoint indexes provided; need relaxation first.
        cls.config_1 = copy.deepcopy(cls.config)
        cls.config_1["common_params"]["is_optimized"] = False
        cls.config_1["common_params"]["wf_name"] = "NEB_test_1"

        # Config 2: The parent structure & two endpoint indexes provided; no need to relax.
        cls.config_2 = copy.deepcopy(cls.config)
        del cls.config_2["fireworks"][0]
        cls.config_2["common_params"]["is_optimized"] = True
        cls.config_2["common_params"]["wf_name"] = "NEB_test_2"

        # Config 3: Two endpoints provided; need to relax two endpoints.
        cls.config_3 = copy.deepcopy(cls.config)
        del cls.config_3["fireworks"][0]
        cls.config_3["common_params"]["is_optimized"] = False
        cls.config_3["common_params"]["wf_name"] = "NEB_test_3"

        # Config 4: Two relaxed endpoints provided; no need to relax two endpoints.
        cls.config_4 = copy.deepcopy(cls.config_3)
        del cls.config_4["fireworks"][0]
        cls.config_4["common_params"]["is_optimized"] = True
        cls.config_4["common_params"]["wf_name"] = "NEB_test_4"

        # Config 5: All images including two endpoints are provided.
        cls.config_5 = copy.deepcopy(cls.config)
        del cls.config_5["fireworks"][0: 2]
        cls.config_5["common_params"]["wf_name"] = "NEB_test_5"

        cls.wf_1 = wf_nudged_elastic_band([parent], parent, cls.config_1)
        cls.wf_2 = wf_nudged_elastic_band([parent], parent, cls.config_2)
        cls.wf_3 = wf_nudged_elastic_band([ep0, ep1], parent, cls.config_3)
        cls.wf_4 = wf_nudged_elastic_band([ep0, ep1], parent, cls.config_4)
        cls.wf_5 = wf_nudged_elastic_band(cls.structures, parent, cls.config_5)

        # Workflow without the config file
        cls.wf_6 = wf_nudged_elastic_band(cls.structures, parent)

    def setUp(self):
        """
        Basic check for scratch directory and launchpad configurations.
        Launchpad will be reset.
        """
        if os.path.exists(self.scratch_dir):
            shutil.rmtree(self.scratch_dir)
        os.makedirs(self.scratch_dir)
        os.chdir(self.scratch_dir)
        try:
            self.lp = LaunchPad.from_file(os.path.join(db_dir, "my_launchpad.yaml"))
            if LAUNCHPAD_RESET:
                self.lp.reset("", require_password=False)
        except:
            raise unittest.SkipTest('Cannot connect to MongoDB! Is the database server running? '
                                    'Are the credentials correct?')

    def tearDown(self):
        if not DEBUG_MODE:
            shutil.rmtree(self.scratch_dir)
            self.lp.reset("", require_password=False)
            # db = self._get_task_database()
            # for coll in db.collection_names():
            #     if coll != "system.indexes":
            #         db[coll].drop()

    def _simulate_vasprun(self, wf):
        """ Run Fake Vasp for testing purpose."""
        test_dir = os.path.abspath(os.path.join(ref_dir, "neb_wf"))
        neb_ref_dirs = {"parent": os.path.join(test_dir, "1"),
                        "ep0": os.path.join(test_dir, "2"),
                        "ep1": os.path.join(test_dir, "3"),
                        "neb1": os.path.join(test_dir, "4"),
                        "neb2": os.path.join(test_dir, "5")}

        return use_fake_vasp(wf, neb_ref_dirs, params_to_check=["ENCUT"])

    def test_wf(self):
        self.wf_1 = self._simulate_vasprun(self.wf_1)
        self.wf_2 = self._simulate_vasprun(self.wf_2)
        self.wf_3 = self._simulate_vasprun(self.wf_3)
        self.wf_4 = self._simulate_vasprun(self.wf_4)
        self.wf_5 = self._simulate_vasprun(self.wf_5)
        self.wf_6 = self._simulate_vasprun(self.wf_6)

        self.lp.add_wf(self.wf_1)
        self.lp.add_wf(self.wf_2)
        self.lp.add_wf(self.wf_3)
        self.lp.add_wf(self.wf_4)
        self.lp.add_wf(self.wf_5)
        self.lp.add_wf(self.wf_6)

        # Use scratch directory as destination directory for testing
        rapidfire(self.lp, fworker=FWorker(env={"run_dest_root": self.scratch_dir}))

        wf = self.lp.get_wf_by_fw_id(1)
        self.assertTrue(all([s == 'COMPLETED' for s in wf.fw_states.values()]))


if __name__ == "__main__":
    t = TestNudgedElasticBandWorkflow()
    t.test_wf()
    t.tearDown()
