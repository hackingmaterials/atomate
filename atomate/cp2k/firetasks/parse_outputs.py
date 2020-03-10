# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import json
import os
import re
from collections import defaultdict
from datetime import datetime

import numpy as np

from monty.json import MontyEncoder, jsanitize
from pydash.objects import has, get

from atomate.vasp.config import DEFUSE_UNSUCCESSFUL
from fireworks import FiretaskBase, FWAction, explicit_serialize
from fireworks.utilities.fw_serializers import DATETIME_HANDLER

from atomate.common.firetasks.glue_tasks import get_calc_loc
from atomate.utils.utils import env_chk, get_meta_from_structure
from atomate.utils.utils import get_logger
from atomate.cp2k.drones import Cp2kDrone
from atomate.cp2k.database import Cp2kCalcDb

__author__ = 'Nicholas Winner'
__email__ = 'nwinner@berkeley.edu'

logger = get_logger(__name__)


@explicit_serialize
class Cp2kToDb(FiretaskBase):
    """
    Enter a cp2k run into the database. Uses current directory unless you
    specify calc_dir or calc_loc.

    Optional params:
        calc_dir (str): path to dir (on current filesystem) that contains cp2k
            output files. Default: use current working directory.
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
            search for the most recent calc_loc with the matching name
        parse_dos (bool): whether to parse the DOS and store in GridFS.
            Defaults to False.
        bandstructure_mode (str): Set to "uniform" for uniform band structure.
            Set to "line" for line mode. If not set, band structure will not
            be parsed.
        additional_fields (dict): dict of additional fields to add
        db_file (str): path to file containing the database credentials.
            Supports env_chk. Default: write data to JSON file.
        fw_spec_field (str): if set, will update the task doc with the contents
            of this key in the fw_spec.
        defuse_unsuccessful (bool): this is a three-way toggle on what to do if
            your job looks OK, but is actually unconverged (either electronic or
            ionic). True -> mark job as COMPLETED, but defuse children.
            False --> do nothing, continue with workflow as normal. "fizzle"
            --> throw an error (mark this job as FIZZLED)
        task_fields_to_push (dict): if set, will update the next Firework/Firetask
            spec using fields from the task document.
            Format: {key : path} -> fw.spec[key] = task_doc[path]
            The path is a full mongo-style path so subdocuments can be referneced
            using dot notation and array keys can be referenced using the index.
            E.g "calcs_reversed.0.output.outar.run_stats"
    """
    optional_params = ["calc_dir", "calc_loc", "parse_dos", "bandstructure_mode",
                       "additional_fields", "db_file", "fw_spec_field", "defuse_unsuccessful",
                       "task_fields_to_push", "parse_chgcar", "parse_aeccar"]

    def run_task(self, fw_spec):
        # get the directory that contains the cp2k dir to parse
        calc_dir = os.getcwd()
        if "calc_dir" in self:
            calc_dir = self["calc_dir"]
        elif self.get("calc_loc"):
            calc_dir = get_calc_loc(self["calc_loc"], fw_spec["calc_locs"])["path"]

        # parse the cp2k directory
        logger.info("PARSING DIRECTORY: {}".format(calc_dir))

        drone = Cp2kDrone(additional_fields=self.get("additional_fields"),
                          parse_dos=self.get("parse_dos", False))

        # assimilate (i.e., parse)
        task_doc = drone.assimilate(calc_dir)

        # Check for additional keys to set based on the fw_spec
        if self.get("fw_spec_field"):
            task_doc.update(fw_spec[self.get("fw_spec_field")])

        # get the database connection
        db_file = env_chk(self.get('db_file'), fw_spec)

        # db insertion or taskdoc dump
        if not db_file:
            with open("task.json", "w") as f:
                f.write(json.dumps(task_doc, default=DATETIME_HANDLER))
        else:
            mmdb = Cp2kCalcDb.from_db_file(db_file, admin=True)
            t_id = mmdb.insert_task(
                task_doc, use_gridfs=self.get("parse_dos", False))
            logger.info("Finished parsing with task_id: {}".format(t_id))

        defuse_children = False
        if task_doc["state"] != "successful":
            defuse_unsuccessful = self.get("defuse_unsuccessful",
                                           DEFUSE_UNSUCCESSFUL)
            if defuse_unsuccessful is True:
                defuse_children = True
            elif defuse_unsuccessful is False:
                pass
            elif defuse_unsuccessful == "fizzle":
                raise RuntimeError(
                    "VaspToDb indicates that job is not successful "
                    "(perhaps your job did not converge within the "
                    "limit of electronic/ionic iterations)!")
            else:
                raise RuntimeError("Unknown option for defuse_unsuccessful: "
                                   "{}".format(defuse_unsuccessful))

        task_fields_to_push = self.get("task_fields_to_push", None)
        update_spec = {}
        if task_fields_to_push:
            if isinstance(task_fields_to_push, dict):
                for key, path_in_task_doc in task_fields_to_push.items():
                    if has(task_doc, path_in_task_doc):
                        update_spec[key] = get(task_doc, path_in_task_doc)
                    else:
                        logger.warn("Could not find {} in task document. Unable to push to next firetask/firework".format(path_in_task_doc))
            else:
                raise RuntimeError("Inappropriate type {} for task_fields_to_push. It must be a "
                                   "dictionary of format: {key: path} where key refers to a field "
                                   "in the spec and path is a full mongo-style path to a "
                                   "field in the task document".format(type(task_fields_to_push)))

        return FWAction(stored_data={"task_id": task_doc.get("task_id", None)},
                        defuse_children=defuse_children, update_spec=update_spec)