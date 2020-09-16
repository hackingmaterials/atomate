# coding: utf-8
"""
Testing for the CalcDb metaclass
"""
import os
import unittest
from collections import defaultdict

import boto3
from maggma.stores import MemoryStore
from moto import mock_s3


__author__ = 'Jimmy Shen <jmmshn@gmail.com>'

from atomate.utils.database import CalcDb
from atomate.utils.utils import get_logger

MODULE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(MODULE_DIR, "..", "..", "common", "test_files")

logger = get_logger(__name__)

class TestToDb(CalcDb):
    def build_indexes(self, indexes=None, background=True):
        pass
    def reset(self):
        pass


class DatabaseTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.testdb = TestToDb.from_db_file(db_dir + "/db_aws.json")

    @classmethod
    def tearDownClass(cls):
        cls.testdb.connection.drop_database(cls.testdb.db_name)

    def test_s3_valid(self):
        with mock_s3():
            conn = boto3.client("s3")
            conn.create_bucket(Bucket="test_bucket")
            index_store = MemoryStore()
            store = self.testdb.get_store('test')
            store.index = index_store
            store.connect()
            store.update([{"fs_id": "mp-1", "data": "111111111110111111"}])
            res = store.query_one({"fs_id" : "mp-1"})
            self.assertEqual(res['fs_id'], "mp-1")
            self.assertEqual(res['data'], "111111111110111111")

    def test_s3_not_valid(self):
        with mock_s3():
            conn = boto3.client("s3")
            conn.create_bucket(Bucket="test_bucket_2")
            index_store = MemoryStore()
            store = self.testdb.get_store('test2')
            store.index = index_store

            with self.assertRaises(Exception):
                store.connect()

