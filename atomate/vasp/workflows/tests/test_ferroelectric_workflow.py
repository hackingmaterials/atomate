import tarfile

from pathlib import Path

from pymatgen import Structure

from fireworks import FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.utils.testing import AtomateTest
from atomate.vasp.powerups import use_fake_vasp, use_potcar_spec
from atomate.vasp.workflows.base.ferroelectric import get_wf_ferroelectric
from atomate.utils.utils import get_a_unique_id

__author__ = "Tess Smidt, Alex Ganose"
__email__ = "blondegeek@gmail.com"

module_dir = Path(__file__).resolve().parent
db_dir = module_dir / "../../../common/test_files"
ref_dir = module_dir / "../../test_files"


class TestFerroelectricWorkflow(AtomateTest):
    def setUp(self, lpad=True):
        super(TestFerroelectricWorkflow, self).setUp(lpad=lpad)

        polar_file = ref_dir / "ferroelectric_wf/BTO_polar_POSCAR"
        non_polar_file = ref_dir / "ferroelectric_wf/BTO_nonpolar_POSCAR"

        self.bto_polar = Structure.from_file(polar_file)
        self.bto_nonpolar = Structure.from_file(non_polar_file)

        self.wfid = "wfid_" + get_a_unique_id()

        self.ferroelectric_config = {
            "vasp_cmd": ">>vasp_cmd<<",
            "db_file": ">>db_file<<",
            "nimages": 2,
            "relax": True,
            "wfid": self.wfid,
            "add_analysis_task": False,
        }

        self.wf = get_wf_ferroelectric(
            self.bto_polar, self.bto_nonpolar, **self.ferroelectric_config
        )

        untar_test_files()

    def tearDown(self):
        super().tearDown()
        delete_folder(ref_dir / "ferroelectric_wf/scratch")

    def _check_run(self, d, mode):

        # Check polar and nonpolar relaxations
        c = d["calcs_reversed"][0]["output"]["structure"]["lattice"]["c"]
        if mode == "_polar_relaxation":
            self.assertAlmostEqual(c, 4.2157, 2)

        if mode == "_nonpolar_relaxation":
            self.assertAlmostEqual(c, 4.0350, 2)

        # Check interpolated structure
        if mode == "_interpolation_1_polarization":
            self.assertAlmostEqual(c, 4.1345, 2)

        # Check that Outcar has needed keys for polarization analysis.
        if "_polarization" in mode and "processing" not in mode:
            outcar = d["calcs_reversed"][0]["output"]["outcar"]

            self.assertIsNotNone(outcar.get("p_ion"))
            self.assertIsNotNone(outcar.get("p_elec"))
            self.assertIsNotNone(outcar.get("p_elec"))
            self.assertIsNotNone(outcar.get("zval_dict"))

        # Check analysis step.
        if mode == "_polarization_post_processing":
            self.assertAlmostEqual(d["polarization_change_norm"], 46.2887527953)

    def test_wf(self):
        self.wf = get_simulated_wf(self.wf)

        # 2 * relax + 3 * polarization = 5
        self.assertEqual(len(self.wf.fws), 5)

        # check VASP parameters on polarization calculation for interpolated
        # structures
        interpolated_polarization_vis = [
            fw.tasks[7]["incar_update"]["lcalcpol"]
            for fw in self.wf.fws
            if "polarization" in fw.name and "interpolation" in fw.name
        ]

        self.assertTrue(all(interpolated_polarization_vis))

        fw_ids = self.lp.add_wf(self.wf)
        rapidfire(self.lp, fworker=FWorker(env={"db_file": db_dir / "db.json"}))

        # Check polar relaxation
        d = self.get_task_collection().find_one(
            {"task_label": "_polar_relaxation"}
        )
        self._check_run(d, "_polar_relaxation")

        # Check nonpolar relaxation
        d = self.get_task_collection().find_one(
            {"task_label": "_nonpolar_relaxation"}
        )
        self._check_run(d, "_nonpolar_relaxation")

        # Check polarization calculations
        all_d = self.get_task_collection().find(
            {"task_label": {"$regex": ".*polarization"}}
        )
        for d in all_d:
            self._check_run(d, d["task_label"])

        # get a fw that can be used to identify the workflow
        fw_id = list(fw_ids.values())[0]

        # check workflow finished without error
        wf = self.lp.get_wf_by_fw_id(fw_id)
        is_completed = [s == "COMPLETED" for s in wf.fw_states.values()]
        self.assertTrue(all(is_completed))


def get_simulated_wf(wf):
    f_dir = ref_dir / "ferroelectric_wf/scratch"
    bto_ref_dirs = {
        "_polar_relaxation": f_dir / "polar_relaxation",
        "_polar_static": f_dir / "polar_static",
        "_polar_polarization": f_dir / "polar_polarization",
        "_nonpolar_relaxation": f_dir / "nonpolar_relaxation",
        "_nonpolar_static": f_dir / "nonpolar_static",
        "_nonpolar_polarization": f_dir / "nonpolar_polarization",
        "_interpolation_1_static": f_dir / "interpolation_1_static",
        "_interpolation_1_polarization": f_dir / "interpolation_1_polarization",
    }

    wf = use_potcar_spec(wf, vasp_to_db_kwargs={"store_volumetric_data": []})
    wf = use_fake_vasp(
        wf, bto_ref_dirs, params_to_check=["ENCUT", "LWAVE"], check_potcar=False
    )

    return wf


def untar_test_files():
    ferro_dir = ref_dir / "ferroelectric_wf"
    t = tarfile.open(ferro_dir / "test_ferroelectric_workflow.gz.tar")
    t.extractall(ferro_dir / "scratch")


def delete_folder(path):
    for sub in path.iterdir():
        if sub.is_dir():
            delete_folder(sub)
        else:
            sub.unlink()
    path.rmdir()
