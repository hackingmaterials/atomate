# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import json
import os
import zlib

import gridfs
import six
from fireworks import FireTaskBase, FWAction
from fireworks.utilities.fw_serializers import DATETIME_HANDLER
from fireworks.utilities.fw_utilities import explicit_serialize
from matgendb.util import get_settings
from monty.json import MontyEncoder
from monty.os.path import zpath
from pymatgen.io.vasp import Vasprun
from pymongo import MongoClient

from matmethods.utils.utils import env_chk
from matmethods.vasp.drones import MMVaspToDbTaskDrone
from matmethods.utils.utils import get_logger

__author__ = 'Anubhav Jain, Kiran Mathew'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov'

logger = get_logger(__name__)


@explicit_serialize
class VaspToDbTask(FireTaskBase):
    """
    Enter a VASP run into the database. By default, the VASP directory is
    assumed to be the current directory.

    Optional params:
        db_file (str): path to file containing the database credentials.
            Supports env_chk. Default: write data to JSON file.
        vasp_dir (str): path to dir (on current filesystem) that contains VASP
            output files. Default: use current working directory.
        vasp_loc (str OR bool): if True will set most recent vasp_loc. If str
            search for the most recent vasp_loc with the matching name
        parse_dos (bool): whether to parse the DOS and store in GridFS.
            Defaults to False.
        bandstructure_mode (str): Set to "uniform" for uniform band structure.
            Set to "line" for line mode. If not set, band structure will not
            be parsed.
        additional_fields (dict): dict of additional fields to add
    """
    optional_params = ["db_file", "vasp_dir", "vasp_loc", "parse_dos",
                       "bandstructure_mode", "additional_fields"]

    def run_task(self, fw_spec):
        # get the directory that contains the VASP dir to parse
        vasp_dir = os.getcwd()
        if "vasp_dir" in self:
            vasp_dir = self["vasp_dir"]
        elif self.get("vasp_loc"):
            if isinstance(self["vasp_loc"], six.string_types):
                for doc in reversed(fw_spec["vasp_locs"]):
                    if doc["name"] == self["vasp_loc_name"]:
                        vasp_dir = doc["path"]
                        break
            else:
                vasp_dir = fw_spec["vasp_locs"][-1]["path"]
        # parse the VASP directory
        logger.info("PARSING DIRECTORY: {}".format(vasp_dir))
        # get the database connection
        db_file = env_chk(self.get('db_file'), fw_spec)

        task_doc = None

        if not db_file:
            drone = MMVaspToDbTaskDrone(simulate_mode=True)
            task_doc = drone.get_task_doc(vasp_dir)
            with open("task.json", "w") as f:
                f.write(json.dumps(task_doc, default=DATETIME_HANDLER))
        else:
            d = get_settings(db_file)
            drone = MMVaspToDbTaskDrone(host=d["host"], port=d["port"],
                                        database=d["database"],
                                        user=d.get("admin_user"),
                                        password=d.get("admin_password"),
                                        collection=d["collection"],
                                        additional_fields=self.get("additional_fields"),
                                        parse_dos=self.get("parse_dos", False), compress_dos=1,
                                        bandstructure_mode=self.get("bandstructure_mode", False),
                                        compress_bs=1)
            t_id, task_doc = drone.assimilate_return_task_doc(vasp_dir)
            logger.info("Finished parsing with task_id: {}".format(t_id))
        return FWAction(stored_data={"task_id": task_doc.get("task_id", None)},
                        defuse_children= (task_doc["state"] != "successful"))
