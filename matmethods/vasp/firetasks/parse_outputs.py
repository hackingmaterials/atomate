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
                                        additional_fields=self.get(
                                            "additional_fields"),
                                        parse_dos=self.get("parse_dos", False),
                                        compress_dos=1)
            t_id, task_doc = drone.assimilate_return_task_doc(vasp_dir)
            logger.info("Finished parsing with task_id: {}".format(t_id))

            if self.get("bandstructure_mode"):
                logger.info("Attempting to parse band structure...")
                if task_doc["state"] != "successful":
                    logger.warn("Skipping band structure insertion; task was not successful!")
                else:
                    # TODO: probably move this into the Drone itself as suggested by Shyam
                    # parse band structure
                    vasprun = Vasprun(
                        zpath(os.path.join(vasp_dir, "vasprun.xml")),
                        parse_eigen=True, parse_projected_eigen=True)
                    bs = vasprun.get_band_structure(
                        line_mode=(self["bandstructure_mode"] == "line"))
                    bs_json = json.dumps(bs.as_dict(), cls=MontyEncoder)
                    bs_compress = zlib.compress(bs_json, 1)

                    # insert band structure data into database
                    conn = MongoClient(d["host"], d["port"])
                    db = conn[d["database"]]
                    if "admin_user" in d:
                        db.authenticate(d["admin_user"], d["admin_password"])
                    tasks = db[d["collection"]]
                    fs = gridfs.GridFS(db, "bandstructure_fs")
                    bs_id = fs.put(bs_compress)
                    tasks.find_one_and_update({"task_id": t_id}, {
                        "$set": {"calcs_reversed.0.bandstructure_fs_id": bs_id,
                                 "calcs_reversed.0.bandstructure_compression":
                                     "zlib"}})
                    logger.info("Finished parsing band structure.")

        return FWAction(stored_data={"task_id": task_doc.get("task_id", None)},
                        defuse_children= (task_doc["state"] != "successful"))


@explicit_serialize
class ToDbTask(FireTaskBase):
    """
    Enter data from a claculation into the database.
    Utilizes a drone to parse the current directory into the DB file to insert.


    Required params:
        drone_cls (AbstractDrone): Class of the drone to conver the data to dict


    Optional params:
        db_file (str): path to file containing the database credentials.
            Supports env_chk. Default: write data to JSON file.
        calc_dir (str): path to dir (on current filesystem) that contains calculation
            output files. Default: use current working directory.
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
            search for the most recent calc_loc with the matching name
        additional_fields (dict): dict of additional fields to add
    """

    required_params = ["drone"]
    optional_params = ["db_file", "calc_dir", "calc_loc", "additional_fields"]

    def run_task(self, fw_spec):
        # get the directory that contains the VASP dir to parse
        calc_dir = os.getcwd()
        if "calc_dir" in self:
            calc_dir = self["calc_dir"]
        elif self.get("calc_loc"):
            if isinstance(self["calc_loc"], six.string_types):
                for doc in reversed(fw_spec["calc_locs"]):
                    if doc["name"] == self["calc_loc_name"]:
                        vasp_dir = doc["path"]
                        break
            else:
                calc_dir = fw_spec["calc_locs"][-1]["path"]
        # parse the Calc directory
        logger.info("PARSING DIRECTORY: {} USING DRONE: {}".format(calc_dir,self['drone'].__class__.__name__))
        # get the database connection
        db_file = env_chk(self.get('db_file'), fw_spec)

        task_doc = None

        if not db_file:
            drone = self['drone'].__class__(simulate_mode=True)
            task_doc = drone.get_task_doc(calc_dir)
            with open("task.json", "w") as f:
                f.write(json.dumps(task_doc, default=DATETIME_HANDLER))
        else:
            d = get_settings(db_file)
            drone = self['drone'].from_db_doc(dbdoc=d,
                                          additional_fields=self.get("additional_fields"),
                                          options = {"parse_dos": self.get("parse_dos", False),
                                                     "compress_dos" : 1})
            t_id, task_doc = drone.assimilate_return_task_doc(calc_dir)
            logger.info("Finished parsing with task_id: {}".format(t_id))


        return FWAction(stored_data={"task_id": task_doc.get("task_id", None)},
                        defuse_children= (task_doc["state"] != "successful"))
