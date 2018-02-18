# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import os
from collections import defaultdict

from pymongo.database import Database

from fireworks import FiretaskBase, Firework, Workflow, explicit_serialize, FWAction

from atomate.utils.utils import env_chk, get_logger, get_mongolike, recursive_get_result, recursive_update, get_database, get_uri

from atomate.utils.testing import AtomateTest

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

MODULE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)))

logger = get_logger(__name__)


@explicit_serialize
class Task1(FiretaskBase):
    def run_task(self, fw_spec):
        print("task1", fw_spec)
        return FWAction(stored_data={"color": "red"})


@explicit_serialize
class Task2(FiretaskBase):
    def run_task(self, fw_spec):
        print("task2", fw_spec)
        return FWAction(stored_data={"color": "yellow"})


class UtilsTests(AtomateTest):

    def setUp(self, lpad=True):

        super(UtilsTests, self).setUp()
        self.fw1 = Firework(Task1())
        self.fw2 = Firework([Task2(), Task2()], parents=self.fw1)
        self.fw3 = Firework(Task1(), parents=self.fw1)


    def test_env_chk(self):
        fw_spec_valid = {"_fw_env": {"hello": "there"}}
        fw_spec_invalid = {}

        self.assertEqual(env_chk("hello", fw_spec_valid), "hello")
        self.assertEqual(env_chk([1, 2, 3], fw_spec_valid), [1, 2, 3])
        self.assertEqual(env_chk(defaultdict(int), fw_spec_valid), defaultdict(int))
        self.assertEqual(env_chk(">>hello<<", fw_spec_valid), "there")
        self.assertRaises(KeyError, env_chk, ">>hello1<<", fw_spec_valid)
        self.assertEqual(env_chk(">>hello1<<", fw_spec_valid, False), None)

        self.assertRaises(KeyError, env_chk, ">>hello<<", fw_spec_invalid)
        self.assertEqual(env_chk(">>hello<<", fw_spec_invalid, False), None)
        self.assertEqual(env_chk(">>hello<<", fw_spec_invalid, False, "fallback"), "fallback")

        self.assertEqual(env_chk(None, fw_spec_valid, False), None)
        self.assertEqual(env_chk(None, fw_spec_valid, False, "fallback"), "fallback")

    def test_get_mongolike(self):
        d = {"a": [{"b": 1}, {"c": {"d": 2}}], "e": {"f": {"g": 3}}, "g": 4, "h": [5, 6]}

        self.assertEqual(get_mongolike(d, "g"), 4)
        self.assertEqual(get_mongolike(d, "e.f.g"), 3)
        self.assertEqual(get_mongolike(d, "a.0.b"), 1)
        self.assertEqual(get_mongolike(d, "a.1.c.d"), 2)
        self.assertEqual(get_mongolike(d, "h.-1"), 6)

    def test_recursive_get_result(self):
        # Basic functionality with dictionary key
        result_dict = {"output": [{}, {"data": [0, 1, 2, 3]}]}
        out = recursive_get_result({"my_data": ">>output.-1.data.2"}, result_dict)
        self.assertEqual(out["my_data"], 2)
        # Basic functionality with attribute
        task1 = Task1()
        out_attr = recursive_get_result({"fw_name": "a>>_fw_name"}, task1)
        self.assertEqual(out_attr["fw_name"], "{{atomate.utils.tests.test_utils.Task1}}")
        # Basic functionality with callable attribute
        out_attr2 = recursive_get_result({"keys": "a>>keys"}, task1)
        self.assertEqual(list(out_attr2["keys"]), [])
        # Testing as_dict functionality
        out_as_dict = recursive_get_result({"fw_name": ">>_fw_name"}, task1)
        self.assertEqual(out_as_dict["fw_name"], "{{atomate.utils.tests.test_utils.Task1}}")

    def test_recursiveupdate(self):
        d = {"a": {"b": 3}, "c": [4]}

        recursive_update(d, {"c": [5]})
        self.assertEqual(d["c"], [5])

        recursive_update(d, {"a": {"b": 5}})
        self.assertEqual(d["a"]["b"], 5)

        recursive_update(d, {"a": {"b": [6]}})
        self.assertEqual(d["a"]["b"], [6])

        recursive_update(d, {"a": {"b": [7]}})
        self.assertEqual(d["a"]["b"], [7])

    def test_get_uri(self):
        self.assertTrue(MODULE_DIR in get_uri(MODULE_DIR))

    def test_get_database(self):
        d = {"host": "localhost", "port": 27017, "database": "atomate_unittest"}

        self.assertTrue(isinstance(get_database(settings=d), Database))

        db = get_database(os.path.join(MODULE_DIR, "db.json"))
        self.assertTrue(isinstance(db, Database))
        self.assertEqual(db.client.address[0], "localhost")
        self.assertEqual(db.name, "atomate_unittest")
