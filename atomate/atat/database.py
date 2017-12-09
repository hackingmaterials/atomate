# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines the database classes.
"""

from atomate.utils.database import CalcDb
from atomate.utils.utils import get_logger

__author__ = 'Matthew Horton'
__credits__ = 'Anubhav Jain, Kiran Mathew, Josua Vieten'
__email__ = 'mkhorton@lbl.gov'
# based on FeffCalcDb

logger = get_logger(__name__)

class McsqsCalcDb(CalcDb):

    def __init__(self, host="localhost", port=27017,
                 database="sqs", collection="tasks",
                 user=None, password=None):
        super(McsqsCalcDb, self).__init__(host, port, database, collection, user, password)

    def build_indexes(self, indexes=None, background=True):
        _indexes = indexes if indexes else ['anonymous_formula']
        for i in _indexes:
            self.collection.create_index(i, background=background)

    def reset(self):
        self.collection.delete_many({})
        self.build_indexes()