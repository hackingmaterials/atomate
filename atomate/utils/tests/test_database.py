# coding: utf-8
"""
Testing for the CalcDb metaclass
"""
import os
import unittest

import boto3
from maggma.stores import MemoryStore
from moto import mock_s3


__author__ = "Jimmy Shen <jmmshn@gmail.com>"

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
            store = self.testdb.get_store("test")
            store.index = index_store
            store.connect()
            store.update([{"fs_id": "mp-1", "data": "111111111110111111"}])
            res = store.query_one({"fs_id": "mp-1"})
            self.assertEqual(res["fs_id"], "mp-1")
            self.assertEqual(res["data"], "111111111110111111")

    def test_s3_not_valid(self):
        with mock_s3():
            conn = boto3.client("s3")
            conn.create_bucket(Bucket="test_bucket_2")
            index_store = MemoryStore()
            store = self.testdb.get_store("test2")
            store.index = index_store

            with self.assertRaises(Exception):
                store.connect()

    def test_maggma_store_names(self):
        with mock_s3():
            conn = boto3.client("s3")
            conn.create_bucket(Bucket="test_bucket")
            index_store = MemoryStore()
            store = self.testdb.get_store("test")
            store.index = index_store
            store.connect()
            self.assertEqual(store.sub_dir, "atomate_test/")
            prefix_db = TestToDb.from_db_file(db_dir + "/db_aws.json")

            prefix_db.maggma_store_prefix = "new_prefix"
            store = prefix_db.get_store("test")
            self.assertEqual(store.sub_dir, "new_prefix_test/")

    def test_uri(self):
        calc_db = TestToDb(
            host_uri="mongodb://localhost:27017",
            database="test_db_name",
            collection="test_collection",
        )
        calc_db.collection.insert_one({"task_id": "mp-1", "data": "12345"})
        self.assertEqual(calc_db.collection.find_one()["data"], "12345")

        with mock_s3():
            conn = boto3.client("s3")
            conn.create_bucket(Bucket="test_bucket")
            uri_db = TestToDb.from_db_file(db_dir + "/db_aws_uri.json")
            store = uri_db.get_store("test")
            self.assertEqual(store.sub_dir, "atomate_test/")

            store.connect()
            store.update([{"fs_id": "mp-1", "data": "111111111110111111"}])
            res = store.query_one({"fs_id": "mp-1"})
            self.assertEqual(res["fs_id"], "mp-1")
            self.assertEqual(res["data"], "111111111110111111")
