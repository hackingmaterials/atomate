# coding: utf-8


"""
This module defines firetasks for parsing and processing the LAMMPS outputfiles to extract useful
information such as the summary of transport properties and insert them into the database.
"""

import os
import json

from fireworks import FiretaskBase, FWAction
from fireworks.utilities.fw_serializers import DATETIME_HANDLER
from fireworks.utilities.fw_utilities import explicit_serialize

from atomate.utils.utils import get_logger
from atomate.common.firetasks.glue_tasks import get_calc_loc
from atomate.utils.utils import env_chk
from atomate.lammps.drones import LammpsDrone
from atomate.lammps.database import LammpsCalcDb

__author__ = 'Kiran Mathew'
__email__ = "kmathew@lbl.gov"

logger = get_logger(__name__)


@explicit_serialize
class LammpsToDB(FiretaskBase):
    """
    Enter a LAMMPS calculation into the database.

    required_params:
        input_filename (str)

    optional_params:
        calc_dir (str): path to dir (on current filesystem) that contains LAMMPS
            output files. Default: use current working directory.
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
            search for the most recent calc_loc with the matching name
        db_file (str): path to file containing the database credentials.
            Supports env_chk. Default: write data to JSON file.
        fw_spec_field:
        data_filename:
        log_filename:
        dump_filenames:
        diffusion_params
        additional_fields:
    """

    required_params = ["input_filename"]

    optional_params = ["calc_dir", "calc_loc", "db_file", "fw_spec_field",
                       "data_filename", "log_filename", "dump_filenames", "diffusion_params",
                       "additional_fields"]

    def run_task(self, fw_spec):

        # get the directory that contains the LAMMPS run parse.
        calc_dir = os.getcwd()
        if "calc_dir" in self:
            calc_dir = self["calc_dir"]
        elif self.get("calc_loc"):
            calc_dir = get_calc_loc(self["calc_loc"], fw_spec["calc_locs"])["path"]

        # parse the directory
        logger.info("PARSING DIRECTORY: {}".format(calc_dir))

        drone = LammpsDrone(additional_fields=self.get("additional_fields"),
                            diffusion_params=self.get("diffusion_params", None))

        task_doc = drone.assimilate(calc_dir, input_filename=self["input_filename"],
                                    log_filename=self.get("log_filename", "log.lammps"),
                                    is_forcefield=self.get("is_forcefield", False),
                                    data_filename=self.get("data_filename", None),
                                    dump_files=self.get("dump_filenames", None))

        # Check for additional keys to set based on the fw_spec
        if self.get("fw_spec_field"):
            task_doc.update(fw_spec[self.get("fw_spec_field")])

        db_file = env_chk(self.get('db_file'), fw_spec)

        # db insertion
        if not db_file:
            with open("task.json", "w") as f:
                f.write(json.dumps(task_doc, default=DATETIME_HANDLER))
        else:
            mmdb = LammpsCalcDb.from_db_file(db_file)
            # insert the task document
            t_id = mmdb.insert(task_doc)
            logger.info("Finished parsing with task_id: {}".format(t_id))

        return FWAction(stored_data={"task_id": task_doc.get("task_id", None)})
