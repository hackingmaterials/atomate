# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines the database classes.
"""

import datetime
import zlib

import gridfs
from pymongo import MongoClient, ASCENDING, DESCENDING

from matmethods.utils.utils import get_logger

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


class MMDb(object):
    """
    Mongodb interface.
    """
    def __init__(self, host="localhost", port=27017, database="vasp", collection="tasks",
                 user=None, password=None):
        self.host = host
        self.database = database
        self.user = user
        self.password = password
        self.port = port
        try:
            connection = MongoClient(self.host, self.port, j=True)
            self.db = connection[self.database]
        except:
            logger.error("Mongodb connection failed")
            raise Exception
        try:
            if self.user:
                self.db.authenticate(self.user, self.password)
        except:
            logger.error("Mongodb authentication failed")
            raise ValueError
        # set counter collection
        if self.db.counter.find({"_id": "taskid"}).count() == 0:
            self.db.counter.insert({"_id": "taskid", "c": 1})
        self.collection = self.db[collection]

    def build(self, indices=None, background=True):
        """
        Build the indices.

        Args:
            indices (list): list of single field indices to be built.
            background (bool): Run in the background or not.

        TODO: make sure that the index building is sensible and check for
            existing indices.
        """
        _indices = indices if indices else ["formula_pretty", "formula_anonymous",
                                            "output.energy", "output.energy_per_atom"]
        self.collection.create_index("task_id", unique=True, background=background)
        # build single field indices
        for i in _indices:
            self.collection.create_index(i, background=background)
        # build compound indices
        for formula in ("formula_pretty", "formula_anonymous"):
            self.collection.create_index([(formula, ASCENDING),
                                          ("output.energy", DESCENDING),
                                          ("completed_at", DESCENDING)],
                                         background=background)
            self.collection.create_index([(formula, ASCENDING),
                                          ("output.energy_per_atom", DESCENDING),
                                          ("completed_at", DESCENDING)],
                                         background=background)

    def insert_gridfs(self, d, root_coll_name="fs", compress=True):
        """
        Insert the given document ot GridFS.

        Args:
            d (dict): the document
            root_coll_name (string): the root collection name
            compress (bool): Whether to compress the data or not

        Returns:
            file id, the type of compression used.
        """
        if compress:
            d = zlib.compress(d, compress)
        fs = gridfs.GridFS(self.db, root_coll_name)
        fs_id = fs.put(d)
        return fs_id, "zlib"

    def insert(self, d, update_duplicates=True):
        """
        Insert the task document ot the database collection.

        Args:
            d (dict): task document
            update_duplicates (bool): whether to update the duplicates

        Returns:
            task_id on successful insertion
        """
        result = self.collection.find_one({"dir_name": d["dir_name"]}, ["dir_name", "task_id"])
        if result is None or update_duplicates:
            d["last_updated"] = datetime.datetime.today()
            if result is None:
                if ("task_id" not in d) or (not d["task_id"]):
                    d["task_id"] = self.db.counter.find_and_modify(query={"_id": "taskid"},
                                                                   update={"$inc": {"c": 1}})["c"]
                logger.info("Inserting {} with taskid = {}".format(d["dir_name"], d["task_id"]))
            elif update_duplicates:
                d["task_id"] = result["task_id"]
                logger.info("Updating {} with taskid = {}".format(d["dir_name"], d["task_id"]))
            self.collection.update({"dir_name": d["dir_name"]}, {"$set": d}, upsert=True)
            return d["task_id"]
        else:
            logger.info("Skipping duplicate {}".format(d["dir_name"]))
            return None
