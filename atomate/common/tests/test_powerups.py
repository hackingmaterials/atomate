import unittest

from fireworks import Firework, ScriptTask, Workflow
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.util.testing import PymatgenTest

from atomate.common.powerups import (
    add_additional_fields_to_taskdocs,
    add_metadata,
    add_priority,
    add_tags,
    powerup_by_kwargs,
    set_queue_adapter,
)
from atomate.utils.utils import get_fws_and_tasks
from atomate.vasp.workflows.base.core import get_wf

__author__ = "Janine George, Guido Petretto"
__email__ = "janine.george@uclouvain.be"


class ModifiedScriptTask(ScriptTask):
    _fw_name = "ModifiedScriptTask"


class TestPowerups(unittest.TestCase):
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
        cls.bsboltz_wf = get_wf(struct_si, "bandstructure_boltztrap.yaml", vis=vis)

    def test_add_priority(self):
        fw1 = Firework([ScriptTask(script=None)], fw_id=-1)
        fw2 = Firework([ScriptTask(script=None)], parents=[fw1], fw_id=-2)
        fw3 = Firework([ScriptTask(script=None)], parents=[fw1], fw_id=-3)

        wf = Workflow([fw1, fw2, fw3])

        wf = add_priority(wf, 4, 8)
        self.assertEqual(wf.id_fw[-1].spec["_priority"], 4)
        self.assertEqual(wf.id_fw[-2].spec["_priority"], 8)
        self.assertEqual(wf.id_fw[-3].spec["_priority"], 8)

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
                    self.assertEqual(t["additional_fields"]["tags"], ["b", "c"])
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
                    self.assertEqual(t["additional_fields"]["tags"], ["foo", "bar"])
                    b_found += 1
                if "VaspToDb" in str(t):
                    self.assertEqual(t["additional_fields"]["tags"], ["foo", "bar"])
                    v_found += 1
        self.assertEqual(b_found, 1)
        self.assertEqual(v_found, 4)

    def test_powerup_by_kwargs(self):
        my_wf = copy_wf(self.bs_wf)
        my_wf = powerup_by_kwargs(
            my_wf,
            [
                {"powerup_name": "add_tags", "kwargs": {"tags_list": ["foo", "bar"]}},
                {
                    "powerup_name": "atomate.common.powerups.add_priority",
                    "kwargs": {"root_priority": 123},
                },
            ],
        )
        self.assertEqual(my_wf.metadata["tags"], ["foo", "bar"])

    def test_set_queue_adapter(self):
        # test fw_name_constraint
        fw1 = Firework([ScriptTask(script=None)], fw_id=-1, name="Firsttask")
        fw2 = Firework(
            [ScriptTask(script=None)], parents=[fw1], fw_id=-2, name="Secondtask"
        )
        fw3 = Firework(
            [ScriptTask(script=None)], parents=[fw1], fw_id=-3, name="Thirdtask"
        )

        wf = Workflow([fw1, fw2, fw3])
        wf = set_queue_adapter(
            wf, {"test": {"test": 1}}, fw_name_constraint="Secondtask"
        )
        self.assertDictEqual(wf.id_fw[-1].spec, {})
        self.assertDictEqual(
            wf.id_fw[-2].spec, {"_queueadapter": {"test": {"test": 1}}}
        )
        self.assertDictEqual(wf.id_fw[-3].spec, {})

        # test task_name_constraint
        fw1 = Firework([ScriptTask(script=None)], fw_id=-1, name="Firsttask")
        fw2 = Firework(
            [ScriptTask(script=None), ModifiedScriptTask(script=None)],
            parents=[fw1],
            fw_id=-2,
            name="Secondtask",
        )
        fw3 = Firework(
            [ScriptTask(script=None)], parents=[fw1], fw_id=-3, name="Thirdtask"
        )

        wf = Workflow([fw1, fw2, fw3])
        wf = set_queue_adapter(
            wf, {"test": {"test": 1}}, task_name_constraint="ModifiedScriptTask"
        )
        self.assertDictEqual(wf.id_fw[-1].spec, {})
        self.assertDictEqual(
            wf.id_fw[-2].spec, {"_queueadapter": {"test": {"test": 1}}}
        )
        self.assertDictEqual(wf.id_fw[-3].spec, {})

    def test_add_additional_fields_to_taskdocs(self):

        my_wf = copy_wf(self.bsboltz_wf)
        meta_dict = {"foo": "bar", "baz": 42}
        my_wf = add_additional_fields_to_taskdocs(my_wf, meta_dict)

        found = 0

        for fw in my_wf.fws:
            for task in fw.tasks:
                if "ToDb" in str(task):
                    for key, val in meta_dict.items():
                        self.assertEqual(task["additional_fields"][key], val)

                    found += 1

        self.assertEqual(found, 5)

    def test_add_metadata(self):
        my_wf = copy_wf(self.bs_wf)
        my_wf.metadata = {"what": "ever"}
        meta_dict = {"foo": "bar", "baz": 42}
        my_wf = add_metadata(my_wf, meta_dict, fw_name_constraint="NonSCFFW")

        self.assertEqual(my_wf.metadata, {"what": "ever", "foo": "bar", "baz": 42})

        for [fw, _] in get_fws_and_tasks(my_wf, fw_name_constraint="NonSCFFW"):
            for key, val in meta_dict.items():
                self.assertEqual(fw.spec[key], val)


def copy_wf(wf):
    return Workflow.from_dict(wf.to_dict())


if __name__ == "__main__":
    unittest.main()
