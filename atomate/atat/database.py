# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

from monty.json import MontyEncoder

"""
This module defines the database classes.
"""

import zlib
import json
from bson import ObjectId

from pymatgen.electronic_structure.bandstructure import BandStructure, BandStructureSymmLine
from pymatgen.electronic_structure.dos import CompleteDos

import gridfs
from pymongo import ASCENDING, DESCENDING

from atomate.utils.database import CalcDb
from fireworks import Firework
from atomate.utils.utils import get_logger

__author__ = 'Kiran Mathew'
__credits__ = 'Anubhav Jain'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


class AtatCalcDb(CalcDb):
    """
    Class to help manage database insertions of Atat drones
    """

    def __init__(self, host="localhost", port=27017, database="SQS", collection="Mcsqs_results", user=None,
                 password=None):
        super(AtatCalcDb, self).__init__(host, port, database, collection, user, password)

    def build_indexes(self, indexes=None, background=True):
        """
        Build the indexes.

        Args:
            indexes (list): list of single field indexes to be built.
            background (bool): Run in the background or not.

        """
        _indices = indexes if indexes else ['disordered',
            'bestsqs',
            'clusters',
            'num_clusters',
            'user_input_settings',
            'objective_function',
            'walltime',
            'atomate_version',
            'mcsqs_version'
            'spacegroup',
            'scaling_matrix',
            'size',
            'last_updated']

        self.collection.create_index("task_id", unique=True, background=background)
        # build single field indexes
        for i in _indices:
            self.collection.create_index(i, background=background)

    def reset(self):
        self.collection.delete_many({})
        self.build_indexes()


