# coding: utf-8


"""
This module defines a base class for derived database classes that store calculation data.
"""

import datetime
from abc import ABCMeta, abstractmethod
from pymongo import MongoClient, ReturnDocument

from monty.json import jsanitize
from monty.serialization import loadfn

from atomate.utils.utils import get_logger
from maggma.stores import S3Store
from maggma.stores import MongoStore

__author__ = 'Kiran Mathew'
__credits__ = 'Anubhav Jain'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


class CalcDb(metaclass=ABCMeta):

    def __init__(self, host, port, database, collection, user, password, **kwargs):
        self.host = host
        self.db_name = database
        self.user = user
        self.password = password
        self.port = int(port)

        # Optional Maggma store for large obj storage
        self._maggma_store_type = None
        self._maggma_login_kwargs = {}
        self.maggma_stores = {}

        try:
            self.connection = MongoClient(host=self.host, port=self.port,
                                          username=self.user,
                                          password=self.password, **kwargs)
            self.db = self.connection[self.db_name]
        except:
            logger.error("Mongodb connection failed")
            raise Exception
        try:
            if self.user:
                self.db.authenticate(self.user, self.password,
                                     source=kwargs.get("authsource", None))
        except:
            logger.error("Mongodb authentication failed")
            raise ValueError
        self.collection = self.db[collection]

        # set counter collection
        if self.db.counter.find({"_id": "taskid"}).count() == 0:
            self.db.counter.insert_one({"_id": "taskid", "c": 0})
            self.build_indexes()

    @abstractmethod
    def build_indexes(self, indexes=None, background=True):
        """
         Build the indexes.

         Args:
             indexes (list): list of single field indexes to be built.
             background (bool): Run in the background or not.
         """
        pass

    def insert(self, d, update_duplicates=True):
        """
        Insert the task document ot the database collection.

        Args:
            d (dict): task document
            update_duplicates (bool): whether to update the duplicates
        """
        result = self.collection.find_one({"dir_name": d["dir_name"]}, ["dir_name", "task_id"])
        if result is None or update_duplicates:
            d["last_updated"] = datetime.datetime.utcnow()
            if result is None:
                if ("task_id" not in d) or (not d["task_id"]):
                    d["task_id"] = self.db.counter.find_one_and_update(
                        {"_id": "taskid"}, {"$inc": {"c": 1}},
                        return_document=ReturnDocument.AFTER)["c"]
                logger.info("Inserting {} with taskid = {}".format(d["dir_name"], d["task_id"]))
            elif update_duplicates:
                d["task_id"] = result["task_id"]
                logger.info("Updating {} with taskid = {}".format(d["dir_name"], d["task_id"]))
            d = jsanitize(d, allow_bson=True)
            self.collection.update_one({"dir_name": d["dir_name"]},
                                       {"$set": d}, upsert=True)
            return d["task_id"]
        else:
            logger.info("Skipping duplicate {}".format(d["dir_name"]))
            return None

    @abstractmethod
    def reset(self):
        pass

    @classmethod
    def from_db_file(cls, db_file, admin=True):
        """
        Create MMDB from database file. File requires host, port, database,
        collection, and optionally admin_user/readonly_user and
        admin_password/readonly_password

        Args:
            db_file (str): path to the file containing the credentials
            admin (bool): whether to use the admin user

        Returns:
            MMDb object
        """
        creds = loadfn(db_file)

        if admin and "admin_user" not in creds and "readonly_user" in creds:
            raise ValueError("Trying to use admin credentials, "
                             "but no admin credentials are defined. "
                             "Use admin=False if only read_only "
                             "credentials are available.")

        if admin:
            user = creds.get("admin_user")
            password = creds.get("admin_password")
        else:
            user = creds.get("readonly_user")
            password = creds.get("readonly_password")

        kwargs = creds.get("mongoclient_kwargs", {})  # any other MongoClient kwargs can go here ...

        if "authsource" in creds:
            kwargs["authsource"] = creds["authsource"]
        else:
            kwargs["authsource"] = creds["database"]

        calc_db = cls(creds["host"], int(creds.get("port", 27017)), creds["database"], creds["collection"],
                      user, password, **kwargs)

        calc_db._maggma_login_kwargs = creds.get("maggma_login", {})
        if "bucket" in calc_db._maggma_login_kwargs:
            calc_db._maggma_store_type = 's3'
            calc_db.get_obj_store = calc_db._get_s3_store
        ## Implement additional maggma stores here as needed

        return calc_db

    def get_obj_store(self, store_name):
        pass

    def _get_s3_store(self, store_name):
        """
        Add a maggma store to this object for storage of large chunk data
        The maggma store will be stored to self.maggma_store[store_name]

        For aws store, all documents will be stored to the same bucket and the store_name will double as the sub_dir name.

        Args:
            store_name: correspond to the the key within calcs_reversed.0 that will be stored
        """
        index_store_ = MongoStore(
            database=self.db_name,
            collection_name=f"atomate_{store_name}_index",
            host=self.host,
            port=self.port,
            username=self.user,
            password=self.password
        )

        store = S3Store(
            index=index_store_,
            sub_dir=f"atomate_{store_name}",
            key="_id",
            **self._maggma_login_kwargs
        )

        self.maggma_stores[store_name] = store
