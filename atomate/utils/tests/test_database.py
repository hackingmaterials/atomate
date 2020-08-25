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

logger = get_logger(__name__)

class TestToDb(CalcDb):
    def build_indexes(self, indexes=None, background=True):
        pass
    def reset(self):
        pass


class DatabaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.testdb = TestToDb.from_db_file(MODULE_DIR + "/db_aws.json")

    def test_s3_valid(self):
        with mock_s3():
            conn = boto3.client("s3")
            conn.create_bucket(Bucket="test_bucket")
            index_store = MemoryStore()
            self.testdb.get_obj_store('test')
            store = self.testdb.maggma_stores['test']
            store.index = index_store
            store.connect()
            store.update([{"task_id": "mp-1", "data": "111111111110111111"}])
            res = store.query_one({"task_id" : "mp-1"})
            self.assertEqual(res['task_id'], "mp-1")
            self.assertEqual(res['data'], "111111111110111111")

    def test_s3_not_valid(self):
        with mock_s3():
            conn = boto3.client("s3")
            conn.create_bucket(Bucket="test_bucket_1")
            index_store = MemoryStore()
            self.testdb.get_obj_store('test')
            store = self.testdb.maggma_stores['test']
            store.index = index_store

            with self.assertRaises(Exception):
                store.connect()



