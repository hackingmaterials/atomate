# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

from monty.json import MontyEncoder, MontyDecoder

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


class Cp2kCalcDb(CalcDb):
    """
    Class to help manage database insertions of cp2k drones
    """

    def __init__(self, host="localhost", port=27017, database="cp2k", collection="tasks", user=None,
                 password=None, **kwargs):
        super(Cp2kCalcDb, self).__init__(host, port, database, collection, user, password, **kwargs)

    def build_indexes(self, indexes=None, background=True):
        """
        Build the indexes.

        Args:
            indexes (list): list of single field indexes to be built.
            background (bool): Run in the background or not.

        TODO: make sure that the index building is sensible and check for
            existing indexes.
        """
        _indices = indexes if indexes else [
            "formula_pretty", "formula_anonymous",
            "output.energy", "output.energy_per_atom", "dir_name"
        ]
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

    def insert_task(self, task_doc, use_gridfs=False):
        """
        Inserts a task document (e.g., as returned by Drone.assimilate()) into the database.
        Handles putting DOS, band structure and charge density into GridFS as needed.
        During testing, a percentage of runs on some clusters had corrupted AECCAR files when even if everything else about the calculation looked OK.
        So we do a quick check here and only record the AECCARs if they are valid

        Args:
            task_doc: (dict) the task document
            use_gridfs (bool) use gridfs for  bandstructures and DOS
        Returns:
            (int) - task_id of inserted document
        """
        dos = None

        # move dos BS and CHGCAR from doc to gridfs
        if use_gridfs and "calcs_reversed" in task_doc:

            if "dos" in task_doc["calcs_reversed"][0]:  # only store idx=0 (last step)
                dos = json.dumps(task_doc["calcs_reversed"][0]["dos"], cls=MontyEncoder)
                del task_doc["calcs_reversed"][0]["dos"]

        # insert the task document
        t_id = self.insert(task_doc)

        # insert the dos into gridfs and update the task document
        if dos:
            dos_gfs_id, compression_type = self.insert_gridfs(dos, "dos_fs", task_id=t_id)
            self.collection.update_one(
                {"task_id": t_id}, {"$set": {"calcs_reversed.0.dos_compression": compression_type}})
            self.collection.update_one({"task_id": t_id}, {"$set": {"calcs_reversed.0.dos_fs_id": dos_gfs_id}})

        return t_id

    def retrieve_task(self, task_id):
        """
        Retrieves a task document and unpacks the band structure and DOS as dict

        Args:
            task_id: (int) task_id to retrieve

        Returns:
            (dict) complete task document with BS + DOS included

        """
        task_doc = self.collection.find_one({"task_id": task_id})
        calc = task_doc["calcs_reversed"][0]
        if 'dos_fs_id' in calc:
            dos = self.get_dos(task_id)
            calc["dos"] = dos.as_dict()
        return task_doc

    def insert_gridfs(self, d, collection="fs", compress=True, oid=None, task_id=None):
        """
        Insert the given document into GridFS.

        Args:
            d (dict): the document
            collection (string): the GridFS collection name
            compress (bool): Whether to compress the data or not
            oid (ObjectId()): the _id of the file; if specified, it must not already exist in GridFS
            task_id(int or str): the task_id to store into the gridfs metadata
        Returns:
            file id, the type of compression used.
        """
        oid = oid or ObjectId()
        compression_type = None

        if compress:
            d = zlib.compress(d.encode(), compress)
            compression_type = "zlib"

        fs = gridfs.GridFS(self.db, collection)
        if task_id:
            # Putting task id in the metadata subdocument as per mongo specs:
            # https://github.com/mongodb/specifications/blob/master/source/gridfs/gridfs-spec.rst#terms
            fs_id = fs.put(d, _id=oid, metadata={"task_id": task_id, "compression": compression_type})
        else:
            fs_id = fs.put(d, _id=oid, metadata={"compression": compression_type})

        return fs_id, compression_type

    def get_band_structure(self, task_id):
        m_task = self.collection.find_one({"task_id": task_id}, {"calcs_reversed": 1})
        fs_id = m_task['calcs_reversed'][0]['bandstructure_fs_id']
        fs = gridfs.GridFS(self.db, 'bandstructure_fs')
        bs_json = zlib.decompress(fs.get(fs_id).read())
        bs_dict = json.loads(bs_json.decode())
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
        dos_dict = json.loads(dos_json.decode())
        return CompleteDos.from_dict(dos_dict)

    def reset(self):
        self.collection.delete_many({})
        self.db.counter.delete_one({"_id": "taskid"})
        self.db.counter.insert_one({"_id": "taskid", "c": 0})
        self.db.dos_fs.files.delete_many({})
        self.db.dos_fs.chunks.delete_many({})
        self.build_indexes()