# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

import json
import os

from fireworks import explicit_serialize, FiretaskBase, FWAction
from fireworks.utilities.fw_serializers import DATETIME_HANDLER
from atomate.utils.utils import get_calc_loc, env_chk
from atomate.vasp.firetasks.parse_outputs import logger

__author__ = 'Shyam Dwaraknath <shyamd@lbl.gov>, Anubhav Jain <ajain@lbl.gov>'


@explicit_serialize
class ToDbTask(FiretaskBase):
    """
    General task to enter data from a calculation into the database.
    Can use any drone to parse the current directory into the DB file to insert.

    Required params:
        drone (AbstractDrone): Drone to convert the data to dict
        mmdb (MMDb): MMDb object

    Optional params:
        db_file (str): path to file containing the database credentials. Supports env_chk.
            Default: write data to JSON file (None).
        calc_dir (str): path to dir (on current filesystem) that contains calculation output files.
            Default: use current working directory.
        calc_loc (str OR bool): if True will set most recent calc_loc. If str search for the most
            recent calc_loc with the matching name
        options (dict): dict of options to pass into the Drone
        additional_fields (dict): dict of additional fields to add
    """

    required_params = ["drone", "mmdb"]
    optional_params = ["db_file", "calc_dir", "calc_loc", "additional_fields", "options"]

    def run_task(self, fw_spec):
        # get the directory that contains the dir to parse
        calc_dir = os.getcwd()
        if "calc_dir" in self:
            calc_dir = self["calc_dir"]
        elif self.get("calc_loc"):
            calc_dir = get_calc_loc(self["calc_loc"], fw_spec["calc_locs"])["path"]

        # parse the calc directory
        logger.info("PARSING DIRECTORY: {} USING DRONE: {}".format(
            calc_dir, self['drone'].__class__.__name__))
        # get the database connection
        db_file = env_chk(self.get('db_file'), fw_spec)

        drone = self['drone'].__class__()
        task_doc = drone.assimilate(calc_dir)
        if not db_file:
            with open("task.json", "w") as f:
                f.write(json.dumps(task_doc, default=DATETIME_HANDLER))
        else:
            mmdb = self["mmdb"]
            db = mmdb.__class__.from_db_file(db_file)
            # insert the task document
            t_id = db.insert(task_doc)
            logger.info("Finished parsing with task_id: {}".format(t_id))

        return FWAction(stored_data={"task_id": task_doc.get("task_id", None)},
                        defuse_children=(task_doc["state"] != "successful"))
