# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines the database classes.
"""

from atomate.utils.database import MMDb
from atomate.utils.utils import get_logger

__author__ = 'Kiran Mathew'
__credits__ = 'Anubhav Jain'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


class MMLammpsDb(MMDb):

    def __init__(self, host="localhost", port=27017, database="lammps", collection="tasks",
                 user=None, password=None):
        super(MMLammpsDb, self).__init__(host, port, database, collection, user, password)

    def build_indexes(self, indexes=None, background=True):
        pass

    def reset(self):
        pass
