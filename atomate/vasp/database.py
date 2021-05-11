"""
This module defines the database classes.
"""

from typing import Any

from monty.json import MontyEncoder
from pymatgen.io.vasp import Chgcar

import zlib
import json
from bson import ObjectId

from pymatgen.electronic_structure.bandstructure import (
    BandStructure,
    BandStructureSymmLine,
)
from pymatgen.electronic_structure.dos import CompleteDos

import gridfs
from pymongo import ASCENDING, DESCENDING

from atomate.utils.database import CalcDb
from atomate.utils.utils import get_logger
from maggma.stores.aws import S3Store
from monty.dev import deprecated

__author__ = "Kiran Mathew"
__credits__ = "Anubhav Jain"
__email__ = "kmathew@lbl.gov"

logger = get_logger(__name__)
# If we use Maggmastores  we will have to initialize a magmma store for each object typl
OBJ_NAMES = (
    "dos",
    "bandstructure",
    "chgcar",
    "locpot",
    "aeccar0",
    "aeccar1",
    "aeccar2",
    "elfcar",
)


class VaspCalcDb(CalcDb):
    """
    Class to help manage database insertions of Vasp drones
    """

    def __init__(
        self,
        host="localhost",
        port=27017,
        database="vasp",
        collection="tasks",
        user=None,
        password=None,
        **kwargs,
    ):
        super().__init__(host, port, database, collection, user, password, **kwargs)

    def build_indexes(self, indexes=None, background=True):
        """
        Build the indexes.

        Args:
            indexes (list): list of single field indexes to be built.
            background (bool): Run in the background or not.

        TODO: make sure that the index building is sensible and check for
            existing indexes.
        """
        _indices = (
            indexes
            if indexes
            else [
                "formula_pretty",
                "formula_anonymous",
                "output.energy",
                "output.energy_per_atom",
                "dir_name",
            ]
        )
        self.collection.create_index("task_id", unique=True, background=background)
        # build single field indexes
        for i in _indices:
            self.collection.create_index(i, background=background)
        # build compound indexes
        for formula in ("formula_pretty", "formula_anonymous"):
            self.collection.create_index(
                [
                    (formula, ASCENDING),
                    ("output.energy", DESCENDING),
                    ("completed_at", DESCENDING),
                ],
                background=background,
            )
            self.collection.create_index(
                [
                    (formula, ASCENDING),
                    ("output.energy_per_atom", DESCENDING),
                    ("completed_at", DESCENDING),
                ],
                background=background,
            )
        # TODO consider sensible index building for the maggma stores

    def insert_task(self, task_doc, use_gridfs=False):
        """
        Inserts a task document (e.g., as returned by Drone.assimilate()) into the database.
        Handles putting DOS, band structure and charge density into GridFS as needed.
        During testing, a percentage of runs on some clusters had corrupted AECCAR files
        when even if everything else about the calculation looked OK.
        So we do a quick check here and only record the AECCARs if they are valid

        Args:
            task_doc (dict): the task document
            use_gridfs (bool): store the data matching OBJ_NAMES to gridfs.
                    if maggma_store_type is set (ex. "s3") this flag will be ignored
        Returns:
            (int) - task_id of inserted document
        """

        big_data_to_store = {}

        def extract_from_calcs_reversed(obj_key):
            """
            Grab the data from calcs_reversed.0.obj_key and store on gridfs directly or some Maggma store
            Args:
                obj_key: Key of the data in calcs_reversed.0 to store
            """
            calcs_r_data = task_doc["calcs_reversed"][0][obj_key]

            # remove the big object from all calcs_reversed
            # this can catch situations were the drone added the data to more than one calc.
            for i_calcs in range(len(task_doc["calcs_reversed"])):
                del task_doc["calcs_reversed"][i_calcs][obj_key]
            return calcs_r_data

        # drop the data from the task_document and keep them in a separate dictionary (big_data_to_store)
        if (
            self._maggma_store_type is not None or use_gridfs
        ) and "calcs_reversed" in task_doc:
            for data_key in OBJ_NAMES:
                if data_key in task_doc["calcs_reversed"][0]:
                    big_data_to_store[data_key] = extract_from_calcs_reversed(data_key)

        # insert the task document
        t_id = self.insert(task_doc)

        if "calcs_reversed" in task_doc:
            # upload the data to a particular location and store the reference to that location in the task database
            for data_key, data_val in big_data_to_store.items():
                fs_di_, compression_type_ = self.insert_object(
                    use_gridfs=use_gridfs,
                    d=data_val,
                    collection=f"{data_key}_fs",
                    task_id=t_id,
                )
                self.collection.update_one(
                    {"task_id": t_id},
                    {
                        "$set": {
                            f"calcs_reversed.0.{data_key}_compression": compression_type_
                        }
                    },
                )
                self.collection.update_one(
                    {"task_id": t_id},
                    {"$set": {f"calcs_reversed.0.{data_key}_fs_id": fs_di_}},
                )
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
        if "bandstructure_fs_id" in calc:
            bs = self.get_band_structure(task_id)
            calc["bandstructure"] = bs.as_dict()
        if "dos_fs_id" in calc:
            dos = self.get_dos(task_id)
            calc["dos"] = dos.as_dict()
        if "chgcar_fs_id" in calc:
            chgcar = self.get_chgcar(task_id)
            calc["chgcar"] = chgcar
        if "aeccar0_fs_id" in calc:
            aeccar = self.get_aeccar(task_id)
            calc["aeccar0"] = aeccar["aeccar0"]
            calc["aeccar2"] = aeccar["aeccar2"]
        return task_doc

    def insert_object(self, use_gridfs, *args, **kwargs):
        """Insert the object into big object storage, try maggma_store if
            it is availible, if not try storing directly to girdfs.

        Args:
            use_gridfs (bool): Whether to store on gridfs if maggma storage is not availible

        Returns:
            fs_id: The id of the stored object
            compression_type: The compress method of the stored object
        """
        if self._maggma_store_type is not None:
            return self.insert_maggma_store(*args, **kwargs)
        elif use_gridfs:
            return self.insert_gridfs(*args, **kwargs)

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

        # always perform the string conversion when inserting directly to gridfs
        d = json.dumps(d, cls=MontyEncoder)
        if compress:
            d = zlib.compress(d.encode(), compress)
            compression_type = "zlib"

        fs = gridfs.GridFS(self.db, collection)
        m_data = {"compression": compression_type}
        if task_id:
            m_data["task_id"] = task_id
        # Putting task id in the metadata subdocument as per mongo specs:
        # https://github.com/mongodb/specifications/blob/master/source/gridfs/gridfs-spec.rst#terms
        fs_id = fs.put(d, _id=oid, metadata=m_data)

        return fs_id, compression_type

    def insert_maggma_store(
        self, d: Any, collection: str, oid: ObjectId = None, task_id: Any = None
    ):
        """
        Insert the given document into a Maggma store, first check if the store is already

        Args:
            data: the document to be stored
            collection (string): the name prefix for the maggma store
            compress (bool): Whether to compress the data or not
            oid (ObjectId()): the _id of the file; if specified, it must not already exist in GridFS
            task_id(int or str): the task_id to store into the gridfs metadata
        Returns:
            file id, the type of compression used.
        """
        oid = oid or str(ObjectId())
        compression_type = None

        doc = {
            "fs_id": oid,
            "maggma_store_type": self.get_store(collection).__class__.__name__,
            "compression": compression_type,
            "data": d,
        }

        search_keys = [
            "fs_id",
        ]
        if task_id is not None:
            search_keys.append("task_id")
            doc["task_id"] = str(task_id)
        elif isinstance(d, dict) and "task_id" in d:
            search_keys.append("task_id")
            doc["task_id"] = str(d["task_id"])

        # make sure the store is availible
        with self.get_store(collection) as store:
            ping_ = store.index._collection.database.command("ping")
            if ping_.get("ok", 0) != 1.0:
                raise ConnectionError(
                    f"Not connected to the index store of {self.__name__}.maggma_store[{collection}]"
                )
            if isinstance(store, S3Store):
                # TODO find some way to ping the aws service
                # ping_ = self._maggma_stores[collection].s3_bucket._name
                pass

            if store.compress:
                compression_type = "zlib"
                doc["compression"] = "zlib"

            store.update([doc], search_keys)

        return oid, compression_type

    def get_data_from_maggma_or_gridfs(self, task_id, key):
        """
        look for a task, then the object of type key associated with that task
        Returns:
            The data stored on object storage, typically a dictionary
        """
        m_task = self.collection.find_one({"task_id": task_id}, {"calcs_reversed": 1})
        fs_id = m_task["calcs_reversed"][0][f"{key}_fs_id"]
        obj_dict = None
        if self._maggma_store_type is not None:
            with self.get_store(f"{key}_fs") as store:
                obj_dict = store.query_one({"fs_id": fs_id})["data"]

        # if the object cannot be found then try using the grid_fs method
        if obj_dict is not None:
            return obj_dict
        else:
            fs = gridfs.GridFS(self.db, f"{key}_fs")
            bs_json = zlib.decompress(fs.get(fs_id).read())
            obj_dict = json.loads(bs_json.decode())
        return obj_dict

    def get_band_structure(self, task_id):
        """
        Read the BS data into a PMG BandStructure or BandStructureSymmLine object

        Args:
            task_id(int or str): the task_id containing the data
        Returns:
            BandStructure or BandStructureSymmLine
        """
        obj_dict = self.get_data_from_maggma_or_gridfs(task_id, key="bandstructure")
        if obj_dict["@class"] == "BandStructure":
            return BandStructure.from_dict(obj_dict)
        elif obj_dict["@class"] == "BandStructureSymmLine":
            return BandStructureSymmLine.from_dict(obj_dict)
        else:
            raise ValueError(
                "Unknown class for band structure! {}".format(obj_dict["@class"])
            )

    def get_dos(self, task_id):
        """
        Read the DOS data into a PMG DOS object

        Args:
            task_id(int or str): the task_id containing the data
        Returns:
            CompleteDos object
        """
        obj_dict = self.get_data_from_maggma_or_gridfs(task_id, key="dos")
        return CompleteDos.from_dict(obj_dict)

    @deprecated("No longer supported, use get_chgcar instead")
    def get_chgcar_string(self, task_id):
        pass

    def get_chgcar(self, task_id):
        """
        Read the CHGCAR data into a PMG Chgcar object
        Args:
            task_id(int or str): the task_id containing the data
        Returns:
            chgcar: Chgcar object
        """
        obj_dict = self.get_data_from_maggma_or_gridfs(task_id, key="chgcar")
        return Chgcar.from_dict(obj_dict)

    def get_aeccar(self, task_id, check_valid=True):
        """
        Read the AECCAR0 + AECCAR2 grid_fs data into a Chgcar object
        Args:
            task_id(int or str): the task_id containing the gridfs metadata
            check_valid (bool): make sure that the aeccar is positive definite
        Returns:
            {"aeccar0" : Chgcar, "aeccar2" : Chgcar}: dict of Chgcar objects
        """

        obj_dict = self.get_data_from_maggma_or_gridfs(task_id, key="aeccar0")
        aeccar0 = Chgcar.from_dict(obj_dict)
        obj_dict = self.get_data_from_maggma_or_gridfs(task_id, key="aeccar2")
        aeccar2 = Chgcar.from_dict(obj_dict)

        if check_valid and (aeccar0.data["total"] + aeccar2.data["total"]).min() < 0:
            ValueError(f"The AECCAR seems to be corrupted for task_id = {task_id}")

        return {"aeccar0": aeccar0, "aeccar2": aeccar2}

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


def put_file_in_gridfs(
    file_path, db, collection_name=None, compress=False, compression_type=None
):
    """
    Helper function to store a file in gridfs.

    Args:
        file_path (str):path to the files that should be saved.
        db (CalcDb): the interface with the database.
        collection_name (str): optionally modify the name of the collection
            with respect to the one included in the db.
        compress (bool): if True the file will be compressed with zlib.
        compression_type (str): if file is already compressed defines the
            compression type to be stored in the metadata.

    Returns:
        ObjectId: the mongodb id of the file that have been saved.
    """

    with open(file_path, "rb") as f:
        data = f.read()

    if compress:
        data = zlib.compress(data, compress)
        compression_type = "zlib"

    if collection_name is None:
        collection_name = db.collection
    fs = gridfs.GridFS(db.db, collection_name)
    fs_id = fs.put(data, metadata={"compression": compression_type})

    return fs_id
