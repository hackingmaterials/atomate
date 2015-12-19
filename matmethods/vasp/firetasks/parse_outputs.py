import json

from fireworks import FireTaskBase
from fireworks.utilities.fw_utilities import explicit_serialize
from matgendb.creator import VaspToDbTaskDrone
from matmethods.utils.utils import env_chk

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


@explicit_serialize
class VaspToDBTask(FireTaskBase):
    """
    Enter a VASP run into the database. By default, the "vasp_locs" key is used to set

    Required params:
        path_to_db_file (str): path to file containing the database credentials. Supports env_chk.

    Optional params:
        dir_to_parse (str): path to dir (on current filesystem) that contains VASP output files.
        vasp_loc_name (str): if set, will search for the most recent vasp_loc with the matching name
    """

    def __init__(self, parameters=None):
        parameters = parameters if parameters else {}
        self.update(parameters)

    def run_task(self, fw_spec):

        # get the directory that contains the VASP dir to parse
        dir_to_parse = None
        if "dir_to_parse" in self:
            dir_to_parse = self["dir_to_parse"]
        elif "vasp_loc_name" in self:
            for doc in reversed(fw_spec["vasp_locs"]):
                if doc["name"] == self["vasp_loc_name"]:
                    dir_to_parse = doc["path"]
                    break
        else:  # default directory to parse is the most recent one in vasp_loc
            dir_to_parse = fw_spec["vasp_locs"][-1]["path"]

        # parse the VASP directory
        print("PARSING DIRECTORY: {}".format(dir_to_parse))
        drone = VaspToDbTaskDrone(simulate_mode=True)
        task_doc = drone.get_task_doc(dir_to_parse)

        print(task_doc)

        """
        # get the database connection
        db_file = env_chk(self.get('path_to_db_file'), fw_spec)

        with open(db_file) as f:
            db_creds = json.load(f)
            # insert into database
        """