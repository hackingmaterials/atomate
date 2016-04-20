# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import unittest
import os
import getpass
from collections import defaultdict

from matmethods.utils.utils import env_chk, get_ssh_connection, get_logger

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

        self.assertEqual(env_chk(None, fw_spec_valid, False), None)

    #@unittest.skipIf(not os.path.exists(os.path.expanduser("~/.ssh/id_rsa")) and not
    #os.path.exists(os.path.expanduser("~/.ssh/authorized_keys")),
    #                 "no '~/.ssh/id_rsa' private key file paramiko test skipped")
    @unittest.skip("paramiko test skipped")
    def test_remote_filesystem(self):
        username = getpass.getuser()
        host = "localhost"
        ssh = get_ssh_connection(username, host, "~/.ssh/id_rsa")
        if ssh:
            logger.debug("paramiko connection OK")
        else:
            logger.error("paramiko connection ERROR. Make sure that passwordless ssh login is setup "
                     "and your private key is in standard location.  e.g. '~/.ssh/id_rsa'")

