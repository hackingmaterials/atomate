import os
from pathlib import Path


from fireworks.core.fworker import FWorker
from fireworks.core.rocket_launcher import rapidfire
from pymatgen.analysis.structure_matcher import StructureMatcher

from atomate.vasp.powerups import use_fake_vasp, use_potcar_spec
from atomate.vasp.workflows.base.electrode import get_ion_insertion_wf
from atomate.utils.testing import AtomateTest

from pymatgen.core import Structure

__author__ = "Jimmy Shen"
__email__ = "jmmshn@gmail.com"

module_dir = Path(__file__).resolve().parent
db_dir = module_dir / "../../../common/test_files"
ref_dir = module_dir / "../../test_files"
wf_dir = ref_dir / "insertion_wf"

VASP_CMD = None  # for fake VASP
DEBUG_MODE = (
    False  # If true, retains the database and output dirs at the end of the test
)


class TestInsertionWorkflow(AtomateTest):
    def setUp(self):
        super().setUp()
        input_output_dirs = ref_dir / "insertion_wf"
        names = os.walk(input_output_dirs).__next__()[1]
        calc_dirs = {n_: input_output_dirs / n_ for n_ in names}
        base_struct = Structure.from_file(wf_dir / "YPO4-static/inputs/POSCAR")
        sm = StructureMatcher(ltol=0.6, stol=0.6, angle_tol=9)
        # Run the workflow with fake vasp
        wf = get_ion_insertion_wf(
            structure=base_struct,
            structure_matcher=sm,
            working_ion="Mg",
            volumetric_data_type="AECCAR",
            db_file=db_dir / "db.json",
            vasp_powerups=[
                {
                    "powerup_name": "add_modify_incar",
                    "kwargs": {"modify_incar_params": {"incar_update": {"KPAR": 8}}},
                },
                {
                    "powerup_name": "use_fake_vasp",
                    "kwargs": {
                        "ref_dirs": calc_dirs,
                        "check_incar": False,
                        "check_kpoints": False,
                        "check_poscar": False,
                        "check_potcar": False,
                    },
                },
                {"powerup_name": "use_potcar_spec", "kwargs": {}},
            ],
            optimizefw_kwargs={"ediffg": -0.05},
        )

        wf_stop_early = get_ion_insertion_wf(
            structure=base_struct,
            structure_matcher=sm,
            working_ion="Mg",
            volumetric_data_type="AECCAR",
            db_file=db_dir / "db.json",
            max_inserted_atoms=1,
            vasp_powerups=[
                {
                    "powerup_name": "add_modify_incar",
                    "kwargs": {"modify_incar_params": {"incar_update": {"KPAR": 8}}},
                },
                {
                    "powerup_name": "use_fake_vasp",
                    "kwargs": {
                        "ref_dirs": calc_dirs,
                        "check_incar": False,
                        "check_kpoints": False,
                        "check_poscar": False,
                        "check_potcar": False,
                    },
                },
                {"powerup_name": "use_potcar_spec", "kwargs": {}},
            ],
            optimizefw_kwargs={"ediffg": -0.05},
        )

        wf = use_fake_vasp(
            wf,
            calc_dirs,
            check_incar=False,
            check_kpoints=False,
            check_poscar=False,
            check_potcar=False,
        )
        wf = use_potcar_spec(wf)
        self.wf = wf

        wf_stop_early = use_fake_vasp(
            wf_stop_early,
            calc_dirs,
            check_incar=False,
            check_kpoints=False,
            check_poscar=False,
            check_potcar=False,
        )
        wf_stop_early = use_potcar_spec(wf_stop_early)
        self.wf_stop_early = wf_stop_early

    def test_has_inserted(self):
        self.lp.add_wf(self.wf_stop_early)
        rapidfire(
            self.lp,
            fworker=FWorker(
                env={
                    "db_file": os.path.join(db_dir, "db.json"),
                    "vasp_cmd": ["echo", "fake"],
                }
            ),
        )
        formula = self.get_task_collection(coll_name="tasks").distinct("formula_pretty")
        self.assertEqual(set(formula), {"Y2Mg(PO4)2", "YPO4"})

        self.lp.add_wf(self.wf)
        rapidfire(
            self.lp,
            fworker=FWorker(
                env={
                    "db_file": os.path.join(db_dir, "db.json"),
                    "vasp_cmd": ["echo", "fake"],
                }
            ),
        )
        # Check that all of the inserted pretty formulas are present
        formula = self.get_task_collection(coll_name="tasks").distinct("formula_pretty")
        self.assertEqual(
            set(formula), {"Y2Mg(PO4)2", "Y2Mg3(PO4)2", "YMg2PO4", "YMgPO4", "YPO4"}
        )
