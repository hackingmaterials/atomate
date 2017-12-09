# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import json
import os

from fireworks import FiretaskBase, FWAction, explicit_serialize
from fireworks.utilities.fw_serializers import DATETIME_HANDLER

from atomate.utils.utils import env_chk
from atomate.common.firetasks.glue_tasks import get_calc_loc
from atomate.utils.utils import get_logger

from atomate.atat.drones import McsqsDrone
from atomate.atat.database import McsqsCalcDb

__author__ = 'Matthew Horton'
__credit__ = 'Kiran Mathew'
__email__ = 'mkhorton@lbl.gov'
# based on FEFF SpectrumToDbTask

logger = get_logger(__name__)


@explicit_serialize
class McsqsToDbTask(FiretaskBase):
    """
    Parse the output of mcsqs and insert it
    into the database.

    Optional_params:
        calc_dir (str): path that contains mcsqs output files,
        uses current directory by default
        calc_loc (str OR bool): if True will set most recent calc_loc,
        if str search for the most recent calc_loc with the matching name
        db_file (str): path to the db file, will output to json by default
    """

    optional_params = ["calc_dir", "calc_loc", "db_file"]

    def run_task(self, fw_spec):
        calc_dir = os.getcwd()
        if "calc_dir" in self:
            calc_dir = self["calc_dir"]
        elif self.get("calc_loc"):
            calc_dir = get_calc_loc(self["calc_loc"], fw_spec["calc_locs"])["path"]

        logger.info("PARSING DIRECTORY: {}".format(calc_dir))

        db_file = env_chk(self.get("db_file"), fw_spec)

        doc = McsqsDrone().assimilate(calc_dir)

        if not db_file:
            with open("mcsqs_task.json", "w") as f:
                f.write(json.dumps(doc, default=DATETIME_HANDLER))
        else:
            db = McsqsCalcDb.from_db_file(db_file, admin=True)
            db.insert(doc)

        logger.info("Finished parsing mcsqs task")

        return FWAction(stored_data={"bestsqs": doc.get("bestsqs", None)})