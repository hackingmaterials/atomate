# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

from matgendb.util import get_database

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
        for t in self._tasks.find({"output.spacegroup.number": {"$type": 2}}, {"task_id": 1, "output": 1}):
            logger.info("Fixing string spacegroup, tid: {}".format(t["task_id"]))
            sg = int(t["output"]["spacegroup"]["number"])
            self._tasks.update_one({"task_id": t["task_id"]},
                                   {"$set": {"output.spacegroup.number": sg}})
        logger.info("FixTasksBuilder finished.")

    def reset(self):
        logger.warn("Cannot reset FixTasksBuilder!")

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
