import os
import unittest
from pathlib import Path

from fireworks.core.fworker import FWorker
from fireworks.core.rocket_launcher import rapidfire
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure

from atomate.utils.testing import AtomateTest
from atomate.vasp.powerups import add_modify_incar, use_fake_vasp, use_potcar_spec
from atomate.vasp.workflows.base.electrode import get_ion_insertion_wf

try:
    from pymatgen.analysis.defects.utils import ChargeInsertionAnalyzer
except ImportError:
    ChargeInsertionAnalyzer = None


__author__ = "Jimmy Shen"
__email__ = "jmmshn@gmail.com"

module_dir = Path(__file__).resolve().parent
db_dir = module_dir / "../../../common/test_files"
ref_dir = module_dir / "../../test_files"
wf_dir = ref_dir / "insertion_wf"


@unittest.skipIf(
    ChargeInsertionAnalyzer is None, "pymatgen.analysis.defects not installed"
)
class TestInsertionWorkflow(AtomateTest):
    def setUp(self):
        super().setUp()
        names = os.walk(wf_dir).__next__()[1]
        calc_dirs = {n_: wf_dir / n_ for n_ in names}
        base_struct = Structure.from_file(wf_dir / "YPO4-static/inputs/POSCAR")
        sm = StructureMatcher(ltol=0.6, stol=0.6, angle_tol=9)

        # Run the workflow with fake VASP
        wf = get_ion_insertion_wf(
            structure=base_struct,
            structure_matcher=sm,
            working_ion="Mg",
            volumetric_data_type="AECCAR",
            db_file=db_dir / "db.json",
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
        wf = add_modify_incar(wf, modify_incar_params={"incar_update": {"KPAR": 8}})
        wf = use_potcar_spec(wf)
        self.wf = wf

        wf_stop_early = get_ion_insertion_wf(
            structure=base_struct,
            structure_matcher=sm,
            working_ion="Mg",
            volumetric_data_type="AECCAR",
            db_file=db_dir / "db.json",
            max_inserted_atoms=1,
            optimizefw_kwargs={"ediffg": -0.05},
        )

        wf_stop_early = use_fake_vasp(
            wf_stop_early,
            calc_dirs,
            check_incar=False,
            check_kpoints=False,
            check_poscar=False,
            check_potcar=False,
        )
        wf = add_modify_incar(wf, modify_incar_params={"incar_update": {"KPAR": 8}})
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
        formulas = self.get_task_collection(coll_name="tasks").distinct(
            "formula_pretty"
        )
        self.assertEqual(set(formulas), {"YPO4"})

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
        formulas = self.get_task_collection(coll_name="tasks").distinct(
            "formula_pretty"
        )
        self.assertEqual(set(formulas), {"YPO4"})
