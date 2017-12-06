# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines a base class for derived database classes that store calculation data.
"""

import datetime
from abc import ABCMeta, abstractmethod
import six
from pymongo import MongoClient, ReturnDocument

from monty.json import jsanitize
from monty.serialization import loadfn

from atomate.utils.utils import get_logger

__author__ = 'Kiran Mathew'
__credits__ = 'Anubhav Jain'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


class CalcDb(six.with_metaclass(ABCMeta)):

    def __init__(self, host, port, database, collection, user, password,
                counter_label="task_id"):
        self.host = host
        self.db_name = database
        self.user = user
        self.password = password
        self.port = int(port)
        self.counter_label = counter_label
        # this is for historical, backwards-compatible reasons
        # due to "task_id" and "taskid" being used in existing databases
        # (we know this makes no sense, sorry!)
        self.counter_label_no_underscore = counter_label.replace("_", "")
        try:
            self.connection = MongoClient(self.host, self.port)
            self.db = self.connection[self.db_name]
        except:
            logger.error("Mongodb connection failed")
            raise Exception
        try:
            if self.user:
                self.db.authenticate(self.user, self.password)
        except:
            logger.error("Mongodb authentication failed")
            raise ValueError
        self.collection = self.db[collection]

        # set counter collection
        if self.db.counter.find({"_id": self.counter_label_no_underscore}).count() == 0:
            self.db.counter.insert_one({"_id": self.counter_label_no_underscore, "c": 0})
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
        result = self.collection.find_one({"dir_name": d["dir_name"]},
                                          ["dir_name", self.counter_label])
        if result is None or update_duplicates:
            d["last_updated"] = datetime.datetime.today()
            if result is None:
                if (self.counter_label not in d) or (not d[self.counter_label]):
                    d[self.counter_label] = self.db.counter.find_one_and_update(
                        {"_id": self.counter_label_no_underscore}, {"$inc": {"c": 1}},
                        return_document=ReturnDocument.AFTER)["c"]
                logger.info("Inserting {} with {} = {}".format(d["dir_name"],
                                                               self.counter_label_no_underscore,
                                                               d["task_id"]))
            elif update_duplicates:
                d[self.counter_label] = result[self.counter_label]
                logger.info("Updating {} with {} = {}".format(d["dir_name"],
                                                              self.counter_label_no_underscore,
                                                              d["task_id"]))
            d = jsanitize(d, allow_bson=True)
            self.collection.update_one({"dir_name": d["dir_name"]},
                                       {"$set": d}, upsert=True)
            return d[self.counter_label]
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

        if admin:
            user = creds.get("admin_user")
            password = creds.get("admin_password")
        else:
            user = creds.get("readonly_user")
            password = creds.get("readonly_password")

        return cls(creds["host"], int(creds["port"]), creds["database"], creds["collection"],
                   user, password)
