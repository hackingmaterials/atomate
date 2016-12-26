# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

from tqdm import tqdm

from matgendb.util import get_database

from atomate.vasp.builders.base import AbstractBuilder


__author__ = 'Alireza Faghanina <albalu@lbl.gov>'

class TagsBuilder(AbstractBuilder):
    def __init__(self, materials_write, tasks_read, update_all=False):
        """
        Starting with an existing materials collection, adds tags from all tasks (if any)
        in the tasks collection.

        Args:
            materials_write (pymongo.collection): materials collection with write access.
            update_all (bool): if true, updates all docs. If false, only updates docs w/o a stability key.
            tasks_read (pymongo.collection): read-only(for safety) tasks collection.

        """
        self._materials = materials_write
        self.update_all = update_all
        self._tasks = tasks_read

    def run(self):
        print("TagsBuilder starting...")
        self._build_indexes()

        # An incorrect incremental update feature:
        # q = {"tags": {"$exists": True}}
        # if not self.update_all:
        #     q["tags"] = {"$exists": False}

        # The incremental update feature is now disabled until I figure out a good/correct way to implement it (AF)
        q = {}

        # TODO: the logic for the query is incorrect. To build incrementally, you need to keep track of which *tasks* have been processed already, and only update materials for which new tasks are available. Just because a materials has tags already in it, doesn't mean that it doesn't require more tags. Or, simply remove the incremental feature for now to prevent incorrect behavior.
        mats = [m for m in self._materials.find(q, {"_tasksbuilder.all_task_ids": 1, "tags": 1,
                                                    "material_id": 1})]
        pbar = tqdm(mats)
        for m in pbar:
            pbar.set_description("Processing materials_id: {}".format(m['material_id']))
            all_tags = []
            try:
                intid_list = [int(taskid.split("-")[1]) for taskid in m["_tasksbuilder"]["all_task_ids"]]
                tasks = self._tasks.find({"task_id": {"$in": intid_list}, "tags": {"$exists": True}},{"tags": 1})
                for task in tasks:
                    all_tags.extend(task["tags"])
                self._materials.update_one({"material_id": m["material_id"]},
                                               {"$set": {"tags": list(set(all_tags))}})
            except:
                import traceback
                print("<---")
                print("There was an error processing material_id: {}, task_id: {}".format(m["material_id"], taskid))
                traceback.print_exc()
                print("--->")
        print("TagsBuilder finished processing.")

    def reset(self):
        self._materials.update_many({}, {"$unset": {"tags": 1}})
        self._build_indexes()

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
