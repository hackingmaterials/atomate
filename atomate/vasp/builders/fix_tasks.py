# coding: utf-8


from atomate.utils.utils import get_database

from atomate.utils.utils import get_logger
from atomate.vasp.builders.base import AbstractBuilder

logger = get_logger(__name__)

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


class FixTasksBuilder(AbstractBuilder):
    def __init__(self, tasks_write):
        """
        Fix historical problems in the tasks database

        Args:
            tasks_write (pymongo.collection): mongodb collection for tasks (write access needed)
        """
        self._tasks = tasks_write

    def run(self):
        # change spacegroup numbers from string to integer where needed
        logger.info("FixTasksBuilder started.")
        for t in self._tasks.find({"output.spacegroup.number": {"$type": 2}},
                                  {"task_id": 1, "output": 1}):
            logger.info("Fixing string spacegroup, tid: {}".format(t["task_id"]))
            sg = int(t["output"]["spacegroup"]["number"])
            self._tasks.update_one({"task_id": t["task_id"]},
                                   {"$set": {"output.spacegroup.number": sg}})

        # change tags from string to list where needed
        for t in self._tasks.find({"tags": {"$exists": True}, "tags.0": {"$exists": False}}, {"task_id": 1, "tags": 1}):
            logger.info("Fixing tag (converting to list), tid: {}".format(t["task_id"]))
            self._tasks.update_one({"task_id": t["task_id"]},
                                   {"$set": {"tags": [t["tags"]]}})

        # fix old (incorrect) delta volume percent
        for t in self._tasks.find({"analysis.delta_volume_percent": {"$exists": True}, "analysis.delta_volume_as_percent": {"$exists": False}}, {"task_id": 1, "analysis": 1}):
            logger.info("Converting delta_volume_percent to be on a percentage scale, tid: {}".format(t["task_id"]))
            self._tasks.update_one({"task_id": t["task_id"]},
                                   {"$set": {"analysis.delta_volume_as_percent": t["analysis"]["delta_volume_percent"] * 100}})

        # remove old (incorrect) delta volume percent
        for t in self._tasks.find(
                {"analysis.delta_volume_percent": {"$exists": True},
                 "analysis.delta_volume_as_percent": {"$exists": True}},
                {"task_id": 1}):
            logger.info("Removing delta_volume_percent, tid: {}".format(t["task_id"]))
            self._tasks.update_one({"task_id": t["task_id"]},
                                   {"$unset": {"analysis.delta_volume_percent": 1}})


        logger.info("FixTasksBuilder finished.")

    def reset(self):
        logger.warning("Cannot reset FixTasksBuilder!")

    @classmethod
    def from_file(cls, db_file, t="tasks", **kwargs):
        """
        Get a FixTasksBuilder using only a db file.

        Args:
            db_file (str): path to db file
            t (str): name of "tasks" collection
            **kwargs: other params to put into FixTasksBuilder
        """
        db_write = get_database(db_file, admin=True)
        return cls(db_write[t], **kwargs)
