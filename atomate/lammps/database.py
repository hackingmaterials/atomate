# coding: utf-8


"""
This module defines the database classes.
"""

import pymongo

from atomate.utils.database import CalcDb
from atomate.utils.utils import get_logger

__author__ = 'Kiran Mathew'
__credits__ = 'Anubhav Jain'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


class LammpsCalcDb(CalcDb):

    def __init__(self, host="localhost", port=27017, database="lammps", collection="tasks",
                 user=None, password=None, **kwargs):
        super(LammpsCalcDb, self).__init__(host, port, database, collection,
                                           user, password, **kwargs)

    def build_indexes(self, indexes=None, background=True):
        indexes = indexes or []
        self.collection.create_index("task_id", unique=True, background=background)
        self.collection.create_index([("completed_at", pymongo.DESCENDING)], background=background)
        for i in indexes:
            self.collection.create_index(i, background=background)

    def reset(self):
        self.collection.delete_many({})
        self.db.counter.delete_one({"_id": "taskid"})
        self.db.counter.insert_one({"_id": "taskid", "c": 0})
        self.build_indexes()
