# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import unittest
from collections import defaultdict

from matmethods.utils.utils import env_chk, get_logger, get_mongolike

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

logger = get_logger(__name__)


class UtilsTests(unittest.TestCase):

    def test_env_chk(self):
        fw_spec_valid = {"_fw_env": {"hello": "there"}}
        fw_spec_invalid = {}

        self.assertEqual(env_chk("hello", fw_spec_valid), "hello")
        self.assertEqual(env_chk([1, 2, 3], fw_spec_valid), [1, 2, 3])
        self.assertEqual(env_chk(defaultdict(int), fw_spec_valid),
                         defaultdict(int))
        self.assertEqual(env_chk(">>hello<<", fw_spec_valid), "there")
        self.assertRaises(KeyError, env_chk, ">>hello1<<", fw_spec_valid)
        self.assertEqual(env_chk(">>hello1<<", fw_spec_valid, False), None)

        self.assertRaises(KeyError, env_chk, ">>hello<<", fw_spec_invalid)
        self.assertEqual(env_chk(">>hello<<", fw_spec_invalid, False), None)
        self.assertEqual(env_chk(">>hello<<", fw_spec_invalid,
                                 False, "fallback"), "fallback")

        self.assertEqual(env_chk(None, fw_spec_valid, False), None)
        self.assertEqual(env_chk(None, fw_spec_valid, False, "fallback"),
                         "fallback")

    def test_get_mongolike(self):
        d = {"a": [{"b": 1}, {"c": {"d": 2}}], "e": {"f": {"g": 3}}, "g": 4}

        self.assertEqual(get_mongolike(d, "g"), 4)
        self.assertEqual(get_mongolike(d, "e.f.g"), 3)
        self.assertEqual(get_mongolike(d, "a.0.b"), 1)
        self.assertEqual(get_mongolike(d, "a.1.c.d"), 2)