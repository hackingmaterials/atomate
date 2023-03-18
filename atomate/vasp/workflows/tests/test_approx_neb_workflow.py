from pathlib import Path

from fireworks import FWorker
from fireworks.core.rocket_launcher import rapidfire
from pymatgen.core import Structure

from atomate.utils.testing import AtomateTest
from atomate.vasp.workflows.base.approx_neb import get_aneb_wf

__author__ = "Ann Rutt"
__email__ = "acrutt@lbl.gov"

module_dir = Path(__file__).resolve().parent
db_dir = module_dir / "../../../common/test_files"
ref_dir = module_dir / "../../test_files"


class TestApproxNEBWorkflow(AtomateTest):
    def setUp(self):
        super().setUp()

        # get base structure file
        input_host_file = ref_dir / "approx_neb_wf/spinel_MnO2_POSCAR"
        self.host_struct = Structure.from_file(input_host_file)

        # get fake vasp powerup inputs
        f_dir = ref_dir / "approx_neb_wf/"
        ref_dirs = {
            "MnO2 host": f_dir / "host",
            "end point: insert Li 0": f_dir / "ep0",
            "end point: insert Li 1": f_dir / "ep1",
            "image 0+1: 0": f_dir / "im0",
            "image 0+1: 1": f_dir / "im1",
            "image 0+1: 2": f_dir / "im2",
        }
        powerup_dicts = [
            {
                "powerup_name": "atomate.vasp.powerups.use_fake_vasp",
                "kwargs": {
                    "ref_dirs": ref_dirs,
                    "check_incar": False,
                    "check_kpoints": False,
                    "check_poscar": True,
                    "check_potcar": False,
                    "clear_inputs": False,
                },
            },
            {"powerup_name": "atomate.vasp.powerups.use_potcar_spec", "kwargs": {}},
        ]

        # get workflow
        self.wf = get_aneb_wf(
            structure=self.host_struct,
            working_ion="Li",
            insert_coords=[[0, 0, 0.5], [0.5, 0.5, 0]],
            insert_coords_combinations=["0+1"],
            n_images=3,
            powerup_dicts=powerup_dicts,
        )

        self.wf_id = self.wf.fws[0].tasks[-1]["approx_neb_wf_uuid"]

    def test_wf(self):
        # 1 host + 2 end point + 1 pathfinder = 4
        self.assertEqual(len(self.wf.fws), 4)

        fw_ids = self.lp.add_wf(self.wf)
        rapidfire(self.lp, fworker=FWorker(env={"db_file": db_dir / "db.json"}))

        # 3 images fws are added after running the workflow
        run_wf = self.lp.get_wf_by_fw_id(list(fw_ids.values())[0])
        self.assertEqual(len(run_wf.fws), 7)

        # check task docs
        host_tds = list(
            self.get_task_collection().find(
                {
                    "approx_neb._source_wf_uuid": self.wf_id,
                    "approx_neb.calc_type": "host",
                }
            )
        )
        self.assertEqual(len(host_tds), 1)
        ep_tds = list(
            self.get_task_collection().find(
                {
                    "approx_neb._source_wf_uuid": self.wf_id,
                    "approx_neb.calc_type": "end_point",
                }
            )
        )
        self.assertEqual(len(ep_tds), 2)
        im_tds = list(
            self.get_task_collection().find(
                {
                    "approx_neb._source_wf_uuid": self.wf_id,
                    "approx_neb.calc_type": "image",
                }
            )
        )
        self.assertEqual(len(im_tds), 3)

        # check approx_neb doc
        aneb_doc = self.get_task_collection(coll_name="approx_neb").find_one(
            {"wf_uuid": self.wf_id}
        )
        for k in ["host", "end_points", "pathfinder", "images"]:
            self.assertIn(k, aneb_doc)
        self.assertIn("output", aneb_doc["host"])
        self.assertEqual(len(aneb_doc["end_points"]), 2)
        for e in aneb_doc["end_points"]:
            self.assertIn("output", e)
        self.assertEqual(len(aneb_doc["pathfinder"]), 1)
        self.assertIn("0+1", aneb_doc["pathfinder"])
        self.assertIn("0+1", aneb_doc["images"])
        self.assertEqual(len(aneb_doc["images"]["0+1"]), 3)
        for i in aneb_doc["images"]["0+1"]:
            self.assertIn("output", i)

        # check workflow finished without error
        is_completed = [s == "COMPLETED" for s in run_wf.fw_states.values()]
        self.assertTrue(all(is_completed))
        self.assertEqual(len(is_completed), 7)

        # 3 images fws are added after running the workflow
        run_wf = self.lp.get_wf_by_fw_id(list(fw_ids.values())[0])
        self.assertEqual(len(run_wf.fws), 7)
