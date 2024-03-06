
from pathlib import Path

from atomate.vasp.workflows.base.lattice_dynamics import get_lattice_dynamics_wf
from fireworks import FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.vasp.powerups import use_fake_vasp, use_potcar_spec
from atomate.utils.testing import AtomateTest

from pymatgen.util.testing import PymatgenTest

__author__ = 'Alex Ganose'
__email__ = 'aganose@lbl.gov'

module_dir = Path(__file__).resolve().parent
db_dir = module_dir / "../../../common/test_files"
ref_dir = module_dir / "../../test_files"
wf_dir = ref_dir / "lattice_dynamics_wf"

test_fworker = FWorker(
    env={"db_file": db_dir / "db.json", "shengbte_cmd": "ShengBTE"}
)


class TestLatticeThermalConductivityWorkflow(AtomateTest):

    def setUp(self, lpad=True):
        super().setUp(lpad=lpad)

        self.struct_si = PymatgenTest.get_structure("Si")

        self.wf = get_lattice_dynamics_wf(
            self.struct_si,
            perturbed_structure_kwargs={"rattle_stds": [0.01]},
            calculate_lattice_thermal_conductivity=True
        )

        # change the cutoffs for fit CSLD
        self.wf.fws[-2].tasks[1].update({"cutoffs": [[5., 3., 3.]]})

        # change shengbte mesh density
        self.wf.fws[-1].tasks[1].update(
            {"control_kwargs": {"reciprocal_density": 500}}
        )

    def test_wf(self):
        self.wf = get_simulated_wf(self.wf)
        self.lp.add_wf(self.wf)
        rapidfire(self.lp, fworker=test_fworker)

        self.assertEqual(len(self.wf.fws), 3)

        # check static calculation
        result = self.get_task_collection().find_one(
            {"task_label": "perturbed structure: i=0; rattle_std=0.010"}
        )
        self.check_run(result, mode="structure optimization")

        wf = self.lp.get_wf_by_fw_id(1)
        self.assertTrue(all([s == 'COMPLETED' for s in wf.fw_states.values()]))

    def check_run(self, d, mode):
        return 1
        # if mode not in ["structure optimization", "static dielectric",
        #                 "raman_0_0.005", "raman analysis"]:
        #     raise ValueError("Invalid mode!")
        #
        # if mode not in ["raman analysis"]:
        #     self.assertEqual(d["formula_pretty"], "Si")
        #     self.assertEqual(d["formula_anonymous"], "A")
        #     self.assertEqual(d["nelements"], 1)
        #     self.assertEqual(d["state"], "successful")
        #     self.assertAlmostEqual(d["calcs_reversed"][0]["output"]["structure"]["lattice"]["a"], 3.867, 2)
        #
        # if mode in ["structure optimization"]:
        #     self.assertAlmostEqual(d["output"]["energy"], -10.850, 2)
        #     self.assertAlmostEqual(d["output"]["energy_per_atom"], -5.425, 2)
        #
        # elif mode in ["static dielectric"]:
        #     epsilon = [[13.23245131, -1.98e-06, -1.4e-06],
        #                [-1.98e-06, 13.23245913, 8.38e-06],
        #                [-1.4e-06, 8.38e-06, 13.23245619]]
        #     np.testing.assert_allclose(epsilon, d["output"]["epsilon_static"], rtol=1e-5)
        #
        # elif mode in ["raman_0_0.005"]:
        #     epsilon = [[13.16509632, 0.00850098, 0.00597267],
        #                [0.00850097, 13.25477303, -0.02979572],
        #                [0.00597267, -0.0297953, 13.28883867]]
        #     np.testing.assert_allclose(epsilon, d["output"]["epsilon_static"], rtol=1e-5)


def get_simulated_wf(wf):
    ref_dirs = {"perturbed structure: i=0; rattle_std=0.010": wf_dir / "0"}

    params_to_check = ["ENCUT", "NSW", "EDIFF"]
    wf = use_potcar_spec(wf)
    wf = use_fake_vasp(
        wf, ref_dirs, params_to_check=params_to_check, check_potcar=False
    )

    return wf

