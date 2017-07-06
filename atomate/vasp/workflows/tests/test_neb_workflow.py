# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import ruamel.yaml as yaml
import os
import unittest
import copy

from fireworks.core.fworker import FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.vasp.powerups import use_fake_vasp
from atomate.vasp.workflows.presets.core import wf_nudged_elastic_band
from atomate.utils.testing import AtomateTest

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
class TestNudgedElasticBandWorkflow(AtomateTest):

    def setUp(self):
        """
        1) Basic check for pymatgen configurations.
        2) Setup all test workflow.
        """
        super(TestNudgedElasticBandWorkflow, self).setUp()
        # Structures used for test:
        parent = PymatgenTest.get_structure("Li2O")
        parent.remove_oxidation_states()
        parent.make_supercell(2)
        ep0, ep1 = get_endpoints_from_index(parent, [0, 1])
        neb_dir = [os.path.join(module_dir, "..", "..", "test_files", "neb_wf", "4", "inputs", "{:02}",
                                "POSCAR").format(i) for i in range(5)]
        self.structures = [Structure.from_file(n) for n in neb_dir]

        # Run fake vasp
        test_yaml = os.path.join(module_dir, "../../test_files/neb_wf/config/neb_unittest.yaml")
        with open(test_yaml, 'r') as stream:
            self.config = yaml.safe_load(stream)
            # Use scratch directory as destination directory for testing
            self.config["common_params"]["_fw_env"] = {"run_dest_root": self.scratch_dir}

        # Config 1: The parent structure & two endpoint indexes provided; need relaxation first.
        self.config_1 = copy.deepcopy(self.config)
        self.config_1["common_params"]["is_optimized"] = False
        self.config_1["common_params"]["wf_name"] = "NEB_test_1"

        # Config 2: The parent structure & two endpoint indexes provided; no need to relax.
        self.config_2 = copy.deepcopy(self.config)
        del self.config_2["fireworks"][0]
        self.config_2["common_params"]["is_optimized"] = True
        self.config_2["common_params"]["wf_name"] = "NEB_test_2"

        # Config 3: Two endpoints provided; need to relax two endpoints.
        self.config_3 = copy.deepcopy(self.config)
        del self.config_3["fireworks"][0]
        self.config_3["common_params"]["is_optimized"] = False
        self.config_3["common_params"]["wf_name"] = "NEB_test_3"

        # Config 4: Two relaxed endpoints provided; no need to relax two endpoints.
        self.config_4 = copy.deepcopy(self.config_3)
        del self.config_4["fireworks"][0]
        self.config_4["common_params"]["is_optimized"] = True
        self.config_4["common_params"]["wf_name"] = "NEB_test_4"

        # Config 5: All images including two endpoints are provided.
        self.config_5 = copy.deepcopy(self.config)
        del self.config_5["fireworks"][0: 2]
        self.config_5["common_params"]["wf_name"] = "NEB_test_5"

        self.wf_1 = wf_nudged_elastic_band([parent], parent, self.config_1)
        self.wf_2 = wf_nudged_elastic_band([parent], parent, self.config_2)
        self.wf_3 = wf_nudged_elastic_band([ep0, ep1], parent, self.config_3)
        self.wf_4 = wf_nudged_elastic_band([ep0, ep1], parent, self.config_4)
        self.wf_5 = wf_nudged_elastic_band(self.structures, parent, self.config_5)

        # Workflow without the config file
        self.wf_6 = wf_nudged_elastic_band(self.structures, parent)

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
    unittest.main()
