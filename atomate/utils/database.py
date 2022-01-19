"""
This module defines a base class for derived database classes that store calculation data.
"""

import datetime
from abc import ABCMeta, abstractmethod

from maggma.stores import MongoStore, MongoURIStore, S3Store
from monty.json import jsanitize
from monty.serialization import loadfn
from pymongo import MongoClient, ReturnDocument
from pymongo.uri_parser import parse_uri

from atomate.utils.utils import get_logger

__author__ = "Kiran Mathew"
__credits__ = "Anubhav Jain"
__email__ = "kmathew@lbl.gov"

logger = get_logger(__name__)


class CalcDb(metaclass=ABCMeta):
    def __init__(
        self,
        host: str = None,
        port: int = None,
        database: str = None,
        collection: str = None,
        user: str = None,
        password: str = None,
        host_uri: str = None,
        maggma_store_kwargs: dict = None,
        maggma_store_prefix: str = "atomate",
        **kwargs,
    ):
        """
        Obeject to handle storing calculation data to MongoDB databases.
        The results of calculations will be parsed by a Drone and
        CalcDb is only responsible for putting that data into the database

        Args:
            host: name of MongoDB host
            port: the port number
            database: the name of the MongoDB database
            collection: the collection were the parsed dictionaries will be stored
            user: MongoDB authentication username
            password: MongoDB authentication password
            host_uri: If the uri designation of the mongodb is provided,
                        other authentication information will be ignored
            maggma_store_kwargs: additional kwargs for mongodb login.
                Currently supports:
                    S3 store kwarges:
                        "bucket" : the S3 bucket where the data is stored
                        "s3_profile" : the S3 profile that contains the login information
                                        typically found at ~/.aws/credentials
                        "compress" : Whether compression is used
                        "endpoint_url" : the url used to access the S3 store
            maggma_store_prefix: when using maggma stores, you can set the prefix string.

            **kwargs:
        """
        if maggma_store_kwargs is None:
            maggma_store_kwargs = {}
        self.maggma_store_prefix = maggma_store_prefix
        self.host = host
        self.db_name = database
        self.user = user
        self.password = password
        self.port = int(port) if port is not None else None
        self.host_uri = host_uri

        self._maggma_store_kwargs = (
            maggma_store_kwargs if maggma_store_kwargs is not None else {}
        )

        self._maggma_store_type = None
        if "bucket" in self._maggma_store_kwargs:
            self._maggma_store_type = "s3"
        # Implement additional maggma stores here as needed

        self._maggma_stores = {}

        if host_uri is not None:
            dd_uri = parse_uri(host_uri)
            if dd_uri["database"] is not None:
                self.db_name = dd_uri["database"]
            else:
                self.host_uri = f"{self.host_uri}/{self.db_name}"

            try:
                self.connection = MongoClient(f"{self.host_uri}")
                self.db = self.connection[self.db_name]
            except Exception:
                logger.error("Mongodb connection failed")
                raise Exception
        else:
            try:
                self.connection = MongoClient(
                    host=self.host,
                    port=self.port,
                    username=self.user,
                    password=self.password,
                    **kwargs,
                )
                self.db = self.connection[self.db_name]
            except Exception:
                logger.error("Mongodb connection failed")
                raise Exception

        self.collection = self.db[collection]

        # set counter collection
        if self.db.counter.count_documents({"_id": "taskid"}) == 0:
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

    def insert(self, d, update_duplicates=True):
        """
        Insert the task document to the database collection.

        Args:
            d (dict): task document
            update_duplicates (bool): whether to update the duplicates
        """
        result = self.collection.find_one(
            {"dir_name": d["dir_name"]}, ["dir_name", "task_id"]
        )
        if result is None or update_duplicates:
            d["last_updated"] = datetime.datetime.utcnow()
            if result is None:
                if ("task_id" not in d) or (not d["task_id"]):
                    d["task_id"] = self.db.counter.find_one_and_update(
                        {"_id": "taskid"},
                        {"$inc": {"c": 1}},
                        return_document=ReturnDocument.AFTER,
                    )["c"]
                logger.info(
                    "Inserting {} with taskid = {}".format(d["dir_name"], d["task_id"])
                )
            elif update_duplicates:
                d["task_id"] = result["task_id"]
                logger.info(
                    "Updating {} with taskid = {}".format(d["dir_name"], d["task_id"])
                )
            d = jsanitize(d, allow_bson=True)
            self.collection.update_one(
                {"dir_name": d["dir_name"]}, {"$set": d}, upsert=True
            )
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

        maggma_kwargs = creds.get("maggma_store", {})
        maggma_prefix = creds.get("maggma_store_prefix", "atomate")
        database = creds.get("database", None)

        kwargs = creds.get(
            "mongoclient_kwargs", {}
        )  # any other MongoClient kwargs can go here ...
        if "host_uri" in creds:
            return cls(
                host_uri=creds["host_uri"],
                database=database,
                collection=creds["collection"],
                maggma_store_kwargs=maggma_kwargs,
                maggma_store_prefix=maggma_prefix,
                **kwargs,
            )

        if admin and "admin_user" not in creds and "readonly_user" in creds:
            raise ValueError(
                "Trying to use admin credentials, "
                "but no admin credentials are defined. "
                "Use admin=False if only read_only "
                "credentials are available."
            )

        if admin:
            user = creds.get("admin_user", "")
            password = creds.get("admin_password", "")
        else:
            user = creds.get("readonly_user", "")
            password = creds.get("readonly_password", "")

        if "authsource" in creds:
            kwargs["authsource"] = creds["authsource"]
        else:
            kwargs["authsource"] = creds["database"]

        return cls(
            host=creds["host"],
            port=int(creds.get("port", 27017)),
            database=creds["database"],
            collection=creds["collection"],
            user=user,
            password=password,
            maggma_store_kwargs=maggma_kwargs,
            maggma_store_prefix=maggma_prefix,
            **kwargs,
        )

    def get_store(self, store_name: str):
        """Get the maggma store with a specific name if it exists, if not create it first.

        Args:
            store_name : name of the store desired
        """
        if store_name not in self._maggma_stores:
            if self._maggma_store_type is None:
                logger.warn(
                    "The maggma store was requested but the maggma store type was not set.  Check your DB_FILE"
                )
                return None
            if self._maggma_store_type == "s3":
                self._maggma_stores[store_name] = self._get_s3_store(store_name)
            # Additional stores can be implemented here
            else:
                raise NotImplementedError("Maggma store type not currently supported.")
        return self._maggma_stores[store_name]

    def _get_s3_store(self, store_name):
        """
        Add a maggma store to this object for storage of large chunk data
        The maggma store will be stored to self.maggma_store[store_name]

        For aws store, all documents will be stored to the same bucket
        and the store_name will double as the sub_dir name.

        Args:
            store_name: correspond to the the key within calcs_reversed.0 that will be stored
        """
        if self.host_uri is not None:
            index_store_ = MongoURIStore(
                uri=self.host_uri,
                database=self.db_name,
                collection_name=f"{self.maggma_store_prefix}_{store_name}_index",
                key="fs_id",
            )
        else:
            index_store_ = MongoStore(
                database=self.db_name,
                collection_name=f"{self.maggma_store_prefix}_{store_name}_index",
                host=self.host,
                port=self.port,
                username=self.user,
                password=self.password,
                key="fs_id",
            )

        store = S3Store(
            index=index_store_,
            sub_dir=f"{self.maggma_store_prefix}_{store_name}",
            key="fs_id",
            **self._maggma_store_kwargs,
        )

        return store
