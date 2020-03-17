# coding: utf-8


from tqdm import tqdm

from atomate.vasp.builders.utils import dbid_to_int, dbid_to_str
from atomate.utils.utils import get_database

from atomate.utils.utils import get_logger
from atomate.vasp.builders.tasks_materials import TasksMaterialsBuilder
from atomate.vasp.builders.base import AbstractBuilder

logger = get_logger(__name__)

__author__ = 'Alireza Faghanina <albalu@lbl.gov>, Anubhav Jain <ajain@lbl.gov>'


class TagsBuilder(AbstractBuilder):
    def __init__(self, materials_write, tasks_read, tasks_prefix="t"):
        """
        Starting with an existing materials collection, searches all its component tasks for
        the "tags" and key in the tasks collection and copies them to the materials collection.
        Thus, the "tags" for a material will be the union of all the tags for its component tasks.

        Args:
            materials_write (pymongo.collection): materials collection with write access.
            tasks_read (pymongo.collection): read-only(for safety) tasks collection.
            tasks_prefix (str): the string prefix for tasks, e.g. "t" for a task_id like "t-132"

        """
        self._materials = materials_write
        self._tasks = tasks_read
        self._tasks_prefix = tasks_prefix

    def run(self):
        logger.info("TagsBuilder starting...")
        self._build_indexes()

        logger.info("Initializing list of all new task_ids to process ...")
        previous_task_ids = []
        for m in self._materials.find({"_tagsbuilder": {"$exists": True}},
                                      {"_tagsbuilder.all_task_ids": 1}):
            previous_task_ids.extend(m["_tagsbuilder"]["all_task_ids"])

        previous_task_ids = [dbid_to_int(t) for t in previous_task_ids]

        q = {"tags": {"$exists": True}, "task_id": {"$nin": previous_task_ids},
             "state": "successful"}

        tasks = [t for t in self._tasks.find(q, {"task_id": 1, "tags": 1})]
        pbar = tqdm(tasks)
        for t in pbar:
            try:
                pbar.set_description("Processing task_id: {}".format(t['task_id']))

                # get the corresponding materials id
                m = self._materials.find_one({"_tasksbuilder.all_task_ids":
                                                  dbid_to_str(self._tasks_prefix, t["task_id"])},
                                             {"material_id": 1, "tags": 1,
                                              "_tagsbuilder": 1})
                if m:
                    all_tags = t["tags"]
                    if "tags" in m and m["tags"]:
                        all_tags.extend(m["tags"])

                    all_tasks = [dbid_to_str(self._tasks_prefix, t["task_id"])]
                    if "_tagsbuilder" in m:
                        all_tasks.extend(m["_tagsbuilder"]["all_task_ids"])

                    all_tags = list(set(all_tags))  # filter duplicates
                    self._materials.update_one({"material_id": m["material_id"]},
                                               {"$set": {"tags": all_tags,
                                                         "_tagsbuilder.all_task_ids": all_tasks}})

            except:
                import traceback
                logger.exception("<---")
                logger.exception("There was an error processing task_id: {}".format(t["task_id"]))
                logger.exception(traceback.format_exc())
                logger.exception("--->")
        logger.info("TagsBuilder finished processing.")

    def reset(self):
        logger.info("Resetting TagsBuilder")
        self._materials.update_many({}, {"$unset": {"tags": 1, "_tagsbuilder": 1}})
        self._build_indexes()
        logger.info("Finished resetting TagsBuilder")

    def _build_indexes(self):
        self._materials.create_index("tags")
        self._materials.create_index("_tagsbuilder.all_task_ids")

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
