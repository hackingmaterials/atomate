import os
from pathlib import Path


from fireworks.core.fworker import FWorker
from fireworks.core.rocket_launcher import rapidfire
from pymatgen.analysis.structure_matcher import StructureMatcher

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
            vasp_powerups={
                "add_modify_incar": {
                    "modify_incar_params": {"incar_update": {"KPAR": 8}}
                },
                "use_fake_vasp": {
                    "ref_dirs": calc_dirs,
                    "check_incar": False,
                    "check_kpoints": False,
                    "check_poscar": False,
                    "check_potcar": False,
                },
            },
            optimizefw_kwargs={"ediffg": -0.05},
        )
        self.wf = wf

    def test_(self):
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
        self.assertEquals(
            set(formula), {"Y2Mg(PO4)2", "Y2Mg3(PO4)2", "YMg2PO4", "YMgPO4", "YPO4"}
        )
