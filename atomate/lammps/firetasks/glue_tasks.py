# coding: utf-8


from fireworks import explicit_serialize, FiretaskBase
import os
from atomate.common.firetasks.glue_tasks import get_calc_loc, CopyFiles
from atomate.utils.utils import env_chk, get_logger
from atomate.lammps.database import LammpsCalcDb
from bson import ObjectId

__author__ = 'Kiran Mathew, Eric Sivonxay'
__email__ = 'kmathew@lbl.gov, esivonxay@lbl.gov'


@explicit_serialize
class CopyPackmolOutputs(CopyFiles):
    """
    Copy files from a previous run directory to the current directory.
    Note: must specify either "calc_loc" or "calc_dir" to indicate the directory
        containing the files to copy.

    Optional params:
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
            search for the most recent calc_loc with the matching name.
        calc_dir (str): path to dir that contains VASP output files.
        filesystem (str): remote filesystem. e.g. username@host
        exclude_files (list): list fo filenames to be excluded when copying.
            NOte: by default nothing is excluded.
    """

    optional_params = ["calc_loc", "calc_dir", "filesystem", "exclude_files"]

    def run_task(self, fw_spec):

        calc_loc = get_calc_loc(self["calc_loc"], fw_spec["calc_locs"]) if self.get("calc_loc") else {}
        exclude_files = self.get("exclude_files", [])

        self.setup_copy(self.get("calc_dir", None), filesystem=self.get("filesystem", None),
                        exclude_files=exclude_files, from_path_dict=calc_loc)
        self.copy_files()

@explicit_serialize
class CopyDeepMDModel(CopyFiles):
    """
    Copy the frozen deepmd model necessary to run
    """

    optional_params = ['model_path']

    def run_task(self, fw_spec):

        model_loc, model_name = os.path.split(self['model_path'])

        self.setup_copy(model_loc, files_to_copy=[model_name])
        self.copy_files()

@explicit_serialize
class DeepMDModelFromDB(FiretaskBase):
    """

    """

    required_params = ["query"]

    optional_params = ["db_file"]

    def run_task(self, fw_spec):

        db_file = env_chk(self["db_file"], fw_spec)
        mmdb = LammpsCalcDb.from_db_file(db_file, admin=True)

        query = self["query"]
        if '_id' in query.keys():
            query['_id'] = ObjectId(query['_id'])

        model_doc = mmdb.db.deepmd_models.find_one(query)

        with open('graph.pb', 'wb') as f:
            f.write(model_doc['model'])