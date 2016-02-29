import json
import os

from fireworks import FireTaskBase
from fireworks.utilities.fw_serializers import DATETIME_HANDLER
from fireworks.utilities.fw_utilities import explicit_serialize
from matgendb.creator import VaspToDbTaskDrone
from matgendb.util import get_settings
from matmethods.utils.utils import env_chk

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


@explicit_serialize
class VaspToDBTask(FireTaskBase):
    """
    Enter a VASP run into the database. By default, the VASP directory is assumed to be the current directory.

    Optional params:
        db_file (str): path to file containing the database credentials. Supports env_chk. Default: write data to JSON file.
        vasp_dir (str): path to dir (on current filesystem) that contains VASP output files. Default: use current working directory.
        vasp_loc (str OR bool): if True will set most recent vasp_loc. If str search for the most recent vasp_loc with the matching name
        additional_fields (dict): dict of additional fields to add
    """

    def run_task(self, fw_spec):

        # get the directory that contains the VASP dir to parse
        vasp_dir = os.getcwd()

        if "vasp_dir" in self:
            vasp_dir = self["vasp_dir"]
        elif self.get("vasp_loc"):
            if isinstance(self["vasp_loc"], basestring):
                for doc in reversed(fw_spec["vasp_locs"]):
                    if doc["name"] == self["vasp_loc_name"]:
                        vasp_dir = doc["path"]
                        break
            else:
                vasp_dir = fw_spec["vasp_locs"][-1]["path"]

        # parse the VASP directory
        print("PARSING DIRECTORY: {}".format(vasp_dir))

        # get the database connection
        db_file = env_chk(self.get('db_file'), fw_spec)

        # TODO: Many important options of VaspToDbTaskDrone are not yet supported

        if not db_file:
            drone = VaspToDbTaskDrone(simulate_mode=True)
            task_doc = drone.get_task_doc(vasp_dir)
            with open("task.json", "w") as f:
                f.write(json.dumps(task_doc, default=DATETIME_HANDLER))

        else:
            d = get_settings(db_file)
            drone = VaspToDbTaskDrone(host=d["host"], port=d["port"],  database=d["database"], user=d.get("admin_user"), password=d.get("admin_password"), collection=d["collection"], additional_fields=self.get("additional_fields"))
            t_id = drone.assimilate(vasp_dir)
            print("Finished parsing with task_id: {}".format(t_id))
