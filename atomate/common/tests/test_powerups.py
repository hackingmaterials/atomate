import unittest

from atomate.common.powerups import set_queue_adapter
from fireworks import Firework, ScriptTask, Workflow

__author__ = "Janine George, Guido Petretto"
__email__ = "janine.george@uclouvain.be"


class ModifiedScriptTask(ScriptTask):
    _fw_name = "ModifiedScriptTask"


class TestPowerups(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
