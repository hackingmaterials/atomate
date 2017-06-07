# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

from tqdm import tqdm

from matgendb.util import get_database

from atomate.utils.utils import get_logger
from atomate.vasp.builders.tasks_materials import TasksMaterialsBuilder
from atomate.vasp.builders.base import AbstractBuilder

logger = get_logger(__name__)

__author__ = 'Alireza Faghanina <albalu@lbl.gov>, Anubhav Jain <ajain@lbl.gov>'


class TagsBuilder(AbstractBuilder):
    def __init__(self, materials_write, tasks_read):
        """
        Starting with an existing materials collection, searches all its component tasks for 
        the "tags" and key in the tasks collection and copies them to the materials collection.
        Thus, the "tags" for a material will be the union of all the tags for its component tasks.

        Args:
            materials_write (pymongo.collection): materials collection with write access.
            tasks_read (pymongo.collection): read-only(for safety) tasks collection.

        """
        self._materials = materials_write
        self._tasks = tasks_read

    def run(self):
        logger.info("TagsBuilder starting...")
        self._build_indexes()

        # TODO: @albalu Build incrementally, taking into account which *tasks* have already been processed -computron
        q = {}
        mats = [m for m in self._materials.find(q, {"_tasksbuilder.all_task_ids": 1, "tags": 1,
                                                    "material_id": 1})]
        pbar = tqdm(mats)
        for m in pbar:
            pbar.set_description("Processing materials_id: {}".format(m['material_id']))
            all_tags = []
            try:
                intid_list = [TasksMaterialsBuilder.tid_to_int(tid) for tid in
                              m["_tasksbuilder"]["all_task_ids"]]
                tasks = self._tasks.find({"task_id": {"$in": intid_list},
                                          "tags": {"$exists": True}}, {"tags": 1})
                for task in tasks:
                    all_tags.extend(task["tags"])
                self._materials.update_one({"material_id": m["material_id"]},
                                           {"$set": {"tags": list(set(all_tags))}})
            except:
                import traceback
                logger.exception("<---")
                logger.exception("There was an error processing material_id: {}".format(m["material_id"]))
                logger.exception(traceback.format_exc())
                logger.exception("--->")
        logger.info("TagsBuilder finished processing.")

    def reset(self):
        logger.info("Resetting TagsBuilder")
        self._materials.update_many({}, {"$unset": {"tags": 1}})
        self._build_indexes()
        logger.info("Finished resetting TagsBuilder")

    def _build_indexes(self):
        self._materials.create_index("tags")

    @classmethod
    def from_file(cls, db_file, m="materials", t="tasks", **kwargs):
        """
        Get a TagsCollector using only a db file.

        Args:
            db_file (str): path to db file
            m (str): name of "materials" collection
            **kwargs: other parameters to feed into the builder, e.g. update_all
        """
        db_write = get_database(db_file, admin=True)
        try:
            db_read = get_database(db_file, admin=False)
            db_read.collection_names()  # throw error if auth failed
        except:
            print("Warning: could not get read-only database; using write creds")
            db_read = get_database(db_file, admin=True)
        return cls(db_write[m], db_read[t], **kwargs)
