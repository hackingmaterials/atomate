from pathlib import Path

import ruamel.yaml as yaml
import unittest
import copy

from fireworks.core.fworker import FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.vasp.powerups import use_fake_vasp, use_potcar_spec
from atomate.vasp.workflows.presets.core import wf_nudged_elastic_band
from atomate.utils.testing import AtomateTest

from pymatgen.core import Structure
from pymatgen.util.testing import PymatgenTest

try:
    from pymatgen_diffusion.neb.io import get_endpoints_from_index
    pmgd = True

except ImportError:
    pmgd = False

__author__ = "Hanmei Tang, Iek-Heng Chu, Alex Ganose"
__email__ = "hat003@eng.ucsd.edu, ihchu@eng.ucsd.edu"

module_dir = Path(__file__).resolve().parent
db_dir = module_dir / "../../../common/test_files"
ref_dir = module_dir / "../../test_files"
wf_dir = ref_dir / "neb_wf"


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
        parent = parent.get_sorted_structure()

        ep0, ep1 = get_endpoints_from_index(parent, [0, 1])
        neb_dir = [wf_dir / "4/inputs/{:02}/POSCAR".format(i) for i in range(5)]
        self.structures = [Structure.from_file(n) for n in neb_dir]

        test_yaml = wf_dir / "config/neb_unittest.yaml"
        with open(test_yaml, "r") as stream:
            self.config = yaml.safe_load(stream)

        # Use scratch directory as destination directory for testing
        env = {"run_dest_root": self.scratch_dir}
        self.config["common_params"]["_fw_env"] = env

        # Config 1: The parent structure & two endpoint indexes provided; need
        # relaxation first.
        self.config_1 = copy.deepcopy(self.config)
        self.config_1["common_params"]["is_optimized"] = False
        self.config_1["common_params"]["wf_name"] = "NEB_test_1"

        # Config 2: The parent structure & two endpoint indexes provided; no
        # need to relax.
        self.config_2 = copy.deepcopy(self.config)
        del self.config_2["fireworks"][0]
        self.config_2["common_params"]["is_optimized"] = True
        self.config_2["common_params"]["wf_name"] = "NEB_test_2"

        # Config 3: Two endpoints provided; need to relax two endpoints.
        self.config_3 = copy.deepcopy(self.config)
        del self.config_3["fireworks"][0]
        self.config_3["common_params"]["is_optimized"] = False
        self.config_3["common_params"]["wf_name"] = "NEB_test_3"

        # Config 4: Two relaxed endpoints provided; no need to relax two
        # endpoints.
        self.config_4 = copy.deepcopy(self.config_3)
        del self.config_4["fireworks"][0]
        self.config_4["common_params"]["is_optimized"] = True
        self.config_4["common_params"]["wf_name"] = "NEB_test_4"

        # Config 5: All images including two endpoints are provided.
        self.config_5 = copy.deepcopy(self.config)
        del self.config_5["fireworks"][0:2]
        self.config_5["common_params"]["wf_name"] = "NEB_test_5"

        self.wf_1 = wf_nudged_elastic_band([parent], parent, self.config_1)
        self.wf_2 = wf_nudged_elastic_band([parent], parent, self.config_2)
        self.wf_3 = wf_nudged_elastic_band([ep0, ep1], parent, self.config_3)
        self.wf_4 = wf_nudged_elastic_band([ep0, ep1], parent, self.config_4)
        self.wf_5 = wf_nudged_elastic_band(self.structures, parent, self.config_5)

        # Workflow without the config file
        self.wf_6 = wf_nudged_elastic_band(self.structures, parent)

    def test_wf(self):
        wf_1 = get_simulated_wf(self.wf_1)
        wf_2 = get_simulated_wf(self.wf_2)
        wf_3 = get_simulated_wf(self.wf_3)
        wf_4 = get_simulated_wf(self.wf_4)
        wf_5 = get_simulated_wf(self.wf_5)
        wf_6 = get_simulated_wf(self.wf_6)

        wf_1_ids = self.lp.add_wf(wf_1)
        wf_2_ids = self.lp.add_wf(wf_2)
        wf_3_ids = self.lp.add_wf(wf_3)
        wf_4_ids = self.lp.add_wf(wf_4)
        wf_5_ids = self.lp.add_wf(wf_5)
        wf_6_ids = self.lp.add_wf(wf_6)

        # get fw ids that can be used to identify the workflows from the DB
        fw_wf_1 = list(wf_1_ids.values())[0]
        fw_wf_2 = list(wf_2_ids.values())[0]
        fw_wf_3 = list(wf_3_ids.values())[0]
        fw_wf_4 = list(wf_4_ids.values())[0]
        fw_wf_5 = list(wf_5_ids.values())[0]
        fw_wf_6 = list(wf_6_ids.values())[0]

        fw_ids = [fw_wf_1, fw_wf_2, fw_wf_3, fw_wf_4, fw_wf_5, fw_wf_6]

        # Use scratch directory as destination directory for testing
        fworker = FWorker(env={"run_dest_root": self.scratch_dir})
        rapidfire(self.lp, fworker=fworker)

        for i in fw_ids:
            wf = self.lp.get_wf_by_fw_id(i)
            is_completed = [s == "COMPLETED" for s in wf.fw_states.values()]
            self.assertTrue(all(is_completed))


def get_simulated_wf(wf):
    dirs = {
        "parent": wf_dir / "1",
        "ep0": wf_dir / "2",
        "ep1": wf_dir / "3",
        "neb1": wf_dir / "4",
        "neb2": wf_dir / "5",
    }

    wf = use_potcar_spec(wf)
    wf = use_fake_vasp(wf, dirs, params_to_check=["ENCUT"], check_potcar=False)
    return wf
