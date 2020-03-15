import unittest

from atomate.utils.utils import get_fws_and_tasks
from fireworks import Firework, ScriptTask, Workflow

from atomate.vasp.powerups import (
    add_priority,
    use_custodian,
    add_trackers,
    add_modify_incar,
    add_small_gap_multiply,
    use_scratch_dir,
    remove_custodian,
    add_tags,
    add_wf_metadata,
    add_modify_potcar,
    add_modify_kpoints,
    clean_up_files,
    set_queue_options,
    use_potcar_spec,
)
from atomate.vasp.workflows.base.core import get_wf

from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.util.testing import PymatgenTest


__author__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"


class TestVaspPowerups(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        struct_si = PymatgenTest.get_structure("Si")
        vis = MPRelaxSet(struct_si, force_gamma=True)
        cls.bs_wf = get_wf(
            struct_si,
            "bandstructure.yaml",
            vis=vis,
            common_params={"vasp_cmd": "test_VASP"},
        )
        cls.bsboltz_wf = get_wf(
            struct_si, "bandstructure_boltztrap.yaml", vis=vis
        )

    def test_add_priority(self):
        fw1 = Firework([ScriptTask(script=None)], fw_id=-1)
        fw2 = Firework([ScriptTask(script=None)], parents=[fw1], fw_id=-2)
        fw3 = Firework([ScriptTask(script=None)], parents=[fw1], fw_id=-3)

        wf = Workflow([fw1, fw2, fw3])

        wf = add_priority(wf, 4, 8)
        self.assertEqual(wf.id_fw[-1].spec["_priority"], 4)
        self.assertEqual(wf.id_fw[-2].spec["_priority"], 8)
        self.assertEqual(wf.id_fw[-3].spec["_priority"], 8)

    def test_custodian_powerups(self):
        my_wf = copy_wf(self.bs_wf)
        my_wf = remove_custodian(my_wf)

        for fw in my_wf.fws:
            task_idx = 1 if "structure optimization" in fw.name else 2
            self.assertTrue("RunVaspDirect" in fw.tasks[task_idx]._fw_name)
            self.assertEqual(fw.tasks[task_idx]["vasp_cmd"], "test_VASP")

        my_wf_double_relax = remove_custodian(copy_wf(self.bs_wf))
        my_wf_double_relax = use_custodian(
            my_wf_double_relax,
            fw_name_constraint="structure optimization",
            custodian_params={"job_type": "double_relaxation_run"},
        )

        for fw in my_wf_double_relax.fws:
            if "structure optimization" in fw.name:
                self.assertTrue("RunVaspCustodian" in fw.tasks[1]._fw_name)
                self.assertEqual(
                    fw.tasks[1]["job_type"], "double_relaxation_run"
                )
            else:
                self.assertTrue("RunVaspDirect" in fw.tasks[2]._fw_name)
                self.assertFalse("job_type" in fw.tasks[2])

    def test_modify_incar(self):
        my_wf = add_modify_incar(
            copy_wf(self.bs_wf),
            {"incar_update": {"NCORE": 1}},
            fw_name_constraint="structure optimization",
        )

        for fw in my_wf.fws:
            if "structure optimization" in fw.name:
                self.assertTrue("ModifyIncar" in fw.tasks[1]._fw_name)
                self.assertEqual(fw.tasks[1]["incar_update"], {"NCORE": 1})
            else:
                for t in fw.tasks:
                    self.assertFalse("ModifyIncar" in t["_fw_name"])

    def test_modify_kpoints(self):
        my_wf = add_modify_kpoints(
            copy_wf(self.bs_wf),
            {"kpoints_update": {"kpts": [[3,4,5]]}},
            fw_name_constraint="structure optimization",
        )

        for fw in my_wf.fws:
            if "structure optimization" in fw.name:
                self.assertTrue("ModifyKpoints" in fw.tasks[1]._fw_name)
                self.assertEqual(fw.tasks[1]["kpoints_update"],
                                 {"kpts": [[3,4,5]]})
            else:
                for t in fw.tasks:
                    self.assertFalse("ModifyKpoints" in t["_fw_name"])

    def test_modify_potcar(self):
        my_wf = add_modify_potcar(
            copy_wf(self.bs_wf),
            {"potcar_symbols": {"Si": "Si_alt"}},
            fw_name_constraint="structure optimization",
        )

        for fw in my_wf.fws:
            if "structure optimization" in fw.name:
                self.assertTrue("ModifyPotcar" in fw.tasks[1]._fw_name)
                self.assertEqual(
                    fw.tasks[1]["potcar_symbols"], {"Si": "Si_alt"}
                )
            else:
                for t in fw.tasks:
                    self.assertFalse("ModifyPotcar" in t["_fw_name"])

    def test_add_trackers(self):
        my_wf = add_trackers(copy_wf(self.bs_wf))

        for fw in my_wf.fws:
            self.assertEqual(len(fw.spec["_trackers"]), 2)

    def test_add_small_gap_multiply(self):
        my_wf = copy_wf(self.bs_wf)
        my_wf = add_small_gap_multiply(my_wf, 0.5, 1.5, "static")
        found = False

        for fw in my_wf.fws:
            if "static" in fw.name:
                for t in fw.tasks:
                    if "WriteVasp" in str(t):
                        self.assertEqual(t["small_gap_multiply"], [0.5, 1.5])
                        found = True

        self.assertEqual(found, True)

    def test_use_scratch_dir(self):
        my_wf = copy_wf(self.bs_wf)
        my_wf = use_custodian(my_wf)
        my_wf = use_scratch_dir(my_wf, ">>scratch_dir<<")
        found = 0

        for fw in my_wf.fws:
            for t in fw.tasks:
                if "RunVaspCustodian" in str(t):
                    self.assertEqual(t["scratch_dir"], ">>scratch_dir<<")
                    found += 1

        self.assertEqual(found, 4)

    def test_add_tags(self):
        my_wf = copy_wf(self.bs_wf)
        my_wf.metadata = {"tags": ["a"]}
        my_wf = add_tags(my_wf, ["b", "c"])

        found = 0

        self.assertEqual(my_wf.metadata["tags"], ["a", "b", "c"])
        for fw in my_wf.fws:
            self.assertEqual(fw.spec["tags"], ["b", "c"])
            for t in fw.tasks:
                if "VaspToDb" in str(t):
                    self.assertEqual(
                        t["additional_fields"]["tags"], ["b", "c"]
                    )
                    found += 1
        self.assertEqual(found, 4)

        my_wf = copy_wf(self.bsboltz_wf)
        my_wf = add_tags(my_wf, ["foo", "bar"])

        v_found = 0
        b_found = 0

        self.assertEqual(my_wf.metadata["tags"], ["foo", "bar"])
        for fw in my_wf.fws:
            self.assertEqual(fw.spec["tags"], ["foo", "bar"])
            for t in fw.tasks:
                if "BoltztrapToDb" in str(t):
                    self.assertEqual(
                        t["additional_fields"]["tags"], ["foo", "bar"]
                    )
                    b_found += 1
                if "VaspToDb" in str(t):
                    self.assertEqual(
                        t["additional_fields"]["tags"], ["foo", "bar"]
                    )
                    v_found += 1
        self.assertEqual(b_found, 1)
        self.assertEqual(v_found, 4)

    def test_queue_options(self):
        my_wf = copy_wf(self.bs_wf)
        my_wf = set_queue_options(my_wf, walltime="00:10:00")
        for fw in my_wf.fws:
            self.assertEqual(fw.spec["_queueadapter"]["walltime"], "00:10:00")

        my_wf = copy_wf(self.bs_wf)
        my_wf = set_queue_options(my_wf, qos="flex")
        for fw in my_wf.fws:
            self.assertEqual(fw.spec["_queueadapter"]["qos"], "flex")

        my_wf = copy_wf(self.bs_wf)
        my_wf = set_queue_options(my_wf, time_min="00:02:00")
        for fw in my_wf.fws:
            self.assertEqual(fw.spec["_queueadapter"]["time_min"], "00:02:00")

        my_wf = copy_wf(self.bs_wf)
        my_wf = set_queue_options(
            my_wf, walltime="00:10:00", time_min="00:02:00", qos="flex"
        )
        for fw in my_wf.fws:
            self.assertEqual(fw.spec["_queueadapter"]["qos"], "flex")
            self.assertEqual(fw.spec["_queueadapter"]["time_min"], "00:02:00")
            self.assertEqual(fw.spec["_queueadapter"]["walltime"], "00:10:00")

    def test_add_wf_metadata(self):
        my_wf = copy_wf(self.bs_wf)
        my_wf = add_wf_metadata(my_wf, PymatgenTest.get_structure("Si"))
        self.assertEqual(my_wf.metadata["nelements"], 1)
        self.assertEqual(my_wf.metadata["formula"], "Si2")

    def test_add_clean_up(self):
        my_wf = clean_up_files(self.bs_wf)
        for fw in my_wf.fws:
            fw_names = [t._fw_name for t in fw.tasks]
            # Raises an error if not in list
            clean_idx = fw_names.index(
                "{{atomate.common.firetasks.glue_tasks.DeleteFiles}}"
            )
            self.assertEqual(
                list(fw.tasks[clean_idx].get("files")), ["WAVECAR*"]
            )

    def test_use_potcar_spec(self):
        wf = copy_wf(self.bs_wf)
        wf = use_potcar_spec(wf)

        idx_list = get_fws_and_tasks(wf, task_name_constraint="WriteVasp")
        self.assertTrue(len(idx_list) > 0)

        for idx_fw, idx_t in idx_list:
            task = wf.fws[idx_fw].tasks[idx_t]
            self.assertTrue(task["potcar_spec"])


def copy_wf(wf):
    return Workflow.from_dict(wf.to_dict())
