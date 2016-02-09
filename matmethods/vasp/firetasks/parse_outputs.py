import json
import os

from fireworks import FireTaskBase
from fireworks.utilities.fw_serializers import DATETIME_HANDLER
from fireworks.utilities.fw_utilities import explicit_serialize
from matgendb.creator import VaspToDbTaskDrone
from matmethods.utils.utils import env_chk

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


@explicit_serialize
class VaspToDBTask(FireTaskBase):
    """
    Enter a VASP run into the database. By default, the VASP directory is assumed to be the current directory.

    Optional params:
        path_to_db_file (str): path to file containing the database credentials. Supports env_chk. If not set, the data is written to a JSON file.
        dir_to_parse (str): path to dir (on current filesystem) that contains VASP output files.
        vasp_loc (str OR bool): if True will set most recent vasp_loc. If str search for the most recent vasp_loc with the matching name
    """

    def run_task(self, fw_spec):

        # get the directory that contains the VASP dir to parse
        dir_to_parse = os.getcwd()

        if "dir_to_parse" in self:
            dir_to_parse = self["dir_to_parse"]
        elif self.get("vasp_loc"):
            if isinstance(self["vasp_loc"], basestring):
                for doc in reversed(fw_spec["vasp_locs"]):
                    if doc["name"] == self["vasp_loc_name"]:
                        dir_to_parse = doc["path"]
                        break
            else:
                dir_to_parse = fw_spec["vasp_locs"][-1]["path"]

        # parse the VASP directory
        print("PARSING DIRECTORY: {}".format(dir_to_parse))
        drone = VaspToDbTaskDrone(simulate_mode=True)
        task_doc = drone.get_task_doc(dir_to_parse)

        print(task_doc)

        # TODO: write to JSON if not path_to_db_file

        # get the database connection
        db_file = env_chk(self.get('path_to_db_file'), fw_spec)

        if not db_file:
            with open("task.json", "w") as f:
                f.write(json.dumps(task_doc, default=DATETIME_HANDLER))

        else:
            with open(db_file) as f:
                db_creds = json.load(f)
            # insert into database