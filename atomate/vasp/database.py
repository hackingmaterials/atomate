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
from atomate.utils.utils import get_logger

__author__ = 'Kiran Mathew'
__credits__ = 'Anubhav Jain'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


class VaspCalcDb(CalcDb):
    """
    Class to help manage database insertions of Vasp drones
    """

    def __init__(self, host="localhost", port=27017, database="vasp", collection="tasks", user=None,
                 password=None):
        super(VaspCalcDb, self).__init__(host, port, database, collection, user, password)

    def build_indexes(self, indexes=None, background=True):
        """
        Build the indexes.

        Args:
            indexes (list): list of single field indexes to be built.
            background (bool): Run in the background or not.

        TODO: make sure that the index building is sensible and check for
            existing indexes.
        """
        _indices = indexes if indexes else ["formula_pretty", "formula_anonymous",
                                            "output.energy", "output.energy_per_atom"]
        self.collection.create_index("task_id", unique=True, background=background)
        # build single field indexes
        for i in _indices:
            self.collection.create_index(i, background=background)
        # build compound indexes
        for formula in ("formula_pretty", "formula_anonymous"):
            self.collection.create_index([(formula, ASCENDING),
                                          ("output.energy", DESCENDING),
                                          ("completed_at", DESCENDING)],
                                         background=background)
            self.collection.create_index([(formula, ASCENDING),
                                          ("output.energy_per_atom", DESCENDING),
                                          ("completed_at", DESCENDING)],
                                         background=background)

    def insert_task(self, task_doc, parse_dos=False, parse_bs=False):
        """
        Inserts a task document (e.g., as returned by Drone.assimilate()) into the database.
        Handles putting DOS and band structure into GridFS as needed.

        Args:
            task_doc: (dict) the task document
            parse_dos: (bool) attempt to parse dos in task_doc and insert into Gridfs
            parse_bs: (bool) attempt to parse bandstructure in task_doc and insert into Gridfs

        Returns:
            (int) - task_id of inserted document
        """

        # insert dos into GridFS
        if parse_dos and "calcs_reversed" in task_doc:
            if "dos" in task_doc["calcs_reversed"][0]:  # only store idx=0 DOS
                dos = json.dumps(task_doc["calcs_reversed"][0]["dos"], cls=MontyEncoder)
                gfs_id, compression_type = self.insert_gridfs(dos, "dos_fs")
                task_doc["calcs_reversed"][0]["dos_compression"] = compression_type
                task_doc["calcs_reversed"][0]["dos_fs_id"] = gfs_id
                del task_doc["calcs_reversed"][0]["dos"]

        # insert band structure into GridFS
        if parse_bs and "calcs_reversed" in task_doc:
            if "bandstructure" in task_doc["calcs_reversed"][0]:  # only store idx=0 BS
                bs = json.dumps(task_doc["calcs_reversed"][0]["bandstructure"], cls=MontyEncoder)
                gfs_id, compression_type = self.insert_gridfs(bs, "bandstructure_fs")
                task_doc["calcs_reversed"][0]["bandstructure_compression"] = compression_type
                task_doc["calcs_reversed"][0]["bandstructure_fs_id"] = gfs_id
                del task_doc["calcs_reversed"][0]["bandstructure"]

        # insert the task document and return task_id
        return self.insert(task_doc)

    def insert_gridfs(self, d, collection="fs", compress=True, oid=None):
        """
        Insert the given document into GridFS.

        Args:
            d (dict): the document
            collection (string): the GridFS collection name
            compress (bool): Whether to compress the data or not
            oid (ObjectId()): the _id of the file; if specified, it must not already exist in GridFS

        Returns:
            file id, the type of compression used.
        """
        oid = oid or ObjectId()
        if compress:
            d = zlib.compress(d.encode(), compress)
        fs = gridfs.GridFS(self.db, collection)
        fs_id = fs.put(d, _id=oid)
        return fs_id, "zlib"

    def get_band_structure(self, task_id):
        m_task = self.collection.find_one({"task_id": task_id}, {"calcs_reversed": 1})
        fs_id = m_task['calcs_reversed'][0]['bandstructure_fs_id']
        fs = gridfs.GridFS(self.db, 'bandstructure_fs')
        bs_json = zlib.decompress(fs.get(fs_id).read())
        bs_dict = json.loads(bs_json)
        if bs_dict["@class"] == "BandStructure":
            return BandStructure.from_dict(bs_dict)
        elif bs_dict["@class"] == "BandStructureSymmLine":
            return BandStructureSymmLine.from_dict(bs_dict)
        else:
            raise ValueError("Unknown class for band structure! {}".format(bs_dict["@class"]))

    def get_dos(self, task_id):
        m_task = self.collection.find_one({"task_id": task_id}, {"calcs_reversed": 1})
        fs_id = m_task['calcs_reversed'][0]['dos_fs_id']
        fs = gridfs.GridFS(self.db, 'dos_fs')
        dos_json = zlib.decompress(fs.get(fs_id).read())
        dos_dict = json.loads(dos_json)
        return CompleteDos.from_dict(dos_dict)

    def reset(self):
        self.collection.delete_many({})
        self.db.counter.delete_one({"_id": "taskid"})
        self.db.counter.insert_one({"_id": "taskid", "c": 0})
        self.db.boltztrap.delete_many({})
        self.db.dos_fs.files.delete_many({})
        self.db.dos_fs.chunks.delete_many({})
        self.db.dos_boltztrap_fs.files.delete_many({})
        self.db.dos_boltztrap_fs.chunks.delete_many({})
        self.db.bandstructure_fs.files.delete_many({})
        self.db.bandstructure_fs.chunks.delete_many({})
        self.build_indexes()


# TODO: @albalu, @matk86, @computron - add BoltztrapCalcDB management here -computron, matk86
