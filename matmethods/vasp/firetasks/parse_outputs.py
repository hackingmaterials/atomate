# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import json
import os

from monty.json import MontyEncoder

from fireworks import FireTaskBase, FWAction
from fireworks.utilities.fw_serializers import DATETIME_HANDLER
from fireworks.utilities.fw_utilities import explicit_serialize
from matgendb.util import get_settings
from matmethods.utils.utils import env_chk, get_calc_loc
from matmethods.utils.utils import get_logger
from matmethods.vasp.database import MMDb
from matmethods.vasp.drones import VaspDrone

__author__ = 'Anubhav Jain, Kiran Mathew, Shyam Dwaraknath'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov, shyamd@lbl.gov'

logger = get_logger(__name__)


@explicit_serialize
class VaspToDbTask(FireTaskBase):
    """
    Enter a VASP run into the database. By default, the VASP directory is
    assumed to be the current directory.

    Optional params:
        calc_dir (str): path to dir (on current filesystem) that contains VASP
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
    """
    optional_params = ["calc_dir", "calc_loc", "parse_dos",
                       "bandstructure_mode", "additional_fields", "db_file"]

    def run_task(self, fw_spec):
        # get the directory that contains the VASP dir to parse
        calc_dir = os.getcwd()
        if "calc_dir" in self:
            calc_dir = self["calc_dir"]
        elif self.get("calc_loc"):
            calc_dir = get_calc_loc(self["calc_loc"], fw_spec["calc_locs"])["path"]

        # parse the VASP directory
        logger.info("PARSING DIRECTORY: {}".format(calc_dir))
        # get the database connection
        db_file = env_chk(self.get('db_file'), fw_spec)

        drone = VaspDrone(additional_fields=self.get("additional_fields"),
                          parse_dos=self.get("parse_dos", False), compress_dos=1,
                          bandstructure_mode=self.get("bandstructure_mode", False), compress_bs=1)
        # assimilate
        task_doc = drone.assimilate(calc_dir)
        # db insertion
        if not db_file:
            with open("task.json", "w") as f:
                f.write(json.dumps(task_doc, default=DATETIME_HANDLER))
        else:
            db_config = get_settings(db_file)
            db = MMDb(host=db_config["host"], port=db_config["port"],
                      database=db_config["database"],
                      user=db_config.get("admin_user"), password=db_config.get("admin_password"),
                      collection=db_config["collection"])

            # insert dos into GridFS
            if self.get("parse_dos") and "calcs_reversed" in task_doc:
                for idx, x in enumerate(task_doc["calcs_reversed"]):
                    if "dos" in task_doc["calcs_reversed"][idx]:
                        if idx == 0:  # only store most recent DOS
                            dos = json.dumps(task_doc["calcs_reversed"][idx]["dos"],
                                             cls=MontyEncoder)
                            gfs_id, compression_type = db.insert_gridfs(dos, "dos_fs")
                            task_doc["calcs_reversed"][idx]["dos_compression"] = compression_type
                            task_doc["calcs_reversed"][idx]["dos_fs_id"] = gfs_id
                        del task_doc["calcs_reversed"][idx]["dos"]

            # insert band structure into GridFS
            if self.get("bandstructure_mode") and "calcs_reversed" in task_doc:
                for idx, x in enumerate(task_doc["calcs_reversed"]):
                    if "bandstructure" in task_doc["calcs_reversed"][idx]:
                        if idx == 0:  # only store most recent band structure
                            bs = json.dumps(task_doc["calcs_reversed"][idx]["bandstructure"],
                                            cls=MontyEncoder)
                            gfs_id, compression_type = db.insert_gridfs(bs, "bandstructure_fs")
                            task_doc["calcs_reversed"][idx][
                                "bandstructure_compression"] = compression_type
                            task_doc["calcs_reversed"][idx]["bandstructure_fs_id"] = gfs_id
                        del task_doc["calcs_reversed"][idx]["bandstructure"]

            # insert the task document
            t_id = db.insert(task_doc)

            logger.info("Finished parsing with task_id: {}".format(t_id))
        return FWAction(stored_data={"task_id": task_doc.get("task_id", None)},
                        defuse_children=(task_doc["state"] != "successful"))


@explicit_serialize
class MakeBandgapCut(FireTaskBase):
    """
    General class for making band gap cuts on calculations that result in whether or not to defuse children Fireworks.
    It possibly should be merged with VaspToDBTask...

    """

    optional_params = ['gt','lt',"calc_dir", "calc_loc", "parse_dos",
                       "bandstructure_mode", "additional_fields", "db_file"]

    def run_task(self, fw_spec):
        # get the directory that contains the VASP dir to parse
        calc_dir = os.getcwd()
        if "calc_dir" in self:
            calc_dir = self["calc_dir"]
        elif self.get("calc_loc"):
            calc_dir = get_calc_loc(self["calc_loc"], fw_spec["calc_locs"])["path"]

        drone = VaspDrone(additional_fields=self.get("additional_fields"),
                          parse_dos=self.get("parse_dos", False), compress_dos=1,
                          bandstructure_mode=self.get("bandstructure_mode", False), compress_bs=1)
        # assimilate
        task_doc = drone.assimilate(calc_dir)

        k = task_doc["calcs_reversed"].keys()
        bandgap = task_doc["calcs_reversed"][k[0]]['output']['bandgap']

        if bandgap < self.get('lt',1.0e10) and bandgap > self.get('gt',0.001):
            return FWAction(stored_data={"bandgap": bandgap})
        else:
            return FWAction(stored_data={"bandgap": bandgap},
                            defuse_children=True)