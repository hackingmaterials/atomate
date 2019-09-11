# coding: utf-8


"""
Lammps Drones.
"""

import os
from datetime import datetime

from pymatgen.apps.borg.hive import AbstractDrone
# from pymatgen.io.lammps.output import LammpsLog, LammpsDump, LammpsRun
# from pymatgen.io.lammps.sets import LammpsInputSet

from atomate.utils.utils import get_uri

from atomate.utils.utils import get_logger

__author__ = 'Brandon Wood, Kiran Mathew'
__email__ = 'b.wood@berkeley.edu'

logger = get_logger(__name__)


class LammpsDrone(AbstractDrone):

    __version__ = 0.1

    schema = {
        "root": {"schema", "dir_name", "input", "output", "last_updated", "state", "completed_at"}
    }

    def __init__(self, additional_fields=None, use_full_uri=True, diffusion_params=None):
        """

        Args:
            additional_fields (dict):
            use_full_uri (bool):
            diffusion_params (dict): parameters to the diffusion_analyzer. If specified a summary
                of diffusion statistics will be added.
        """
        self.additional_fields = additional_fields or {}
        self.use_full_uri = use_full_uri
        self.runs = []
        self.diffusion_params = diffusion_params

    def assimilate(self, path, input_filename, log_filename="log.lammps",  is_forcefield=False,
                   data_filename=None, dump_files=None):
        """
        Parses lammps input, data and log files and insert the result into the db.

        Args:
            path (str): path to the run folder
            input_filename (str): just the name of the input file
            log_filename (str): lammps log file name
            is_forcefield (bool): whether or not ot parse forcefield info
            data_filename (str): name of the data file
            dump_files ([str]): list of dump file names

        Returns:
            dict
        """

        input_file = os.path.join(path, input_filename)
        data_file = os.path.join(path, data_filename) if data_filename else None
        log_file = os.path.join(path, log_filename)
        dump_files = dump_files or []
        dump_files = [dump_files] if isinstance(dump_files, str) else dump_files

        # input set
        lmps_input = LammpsInputSet.from_file("lammps", input_file, {}, data_file, data_filename)

        # dumps
        dumps = []
        if dump_files:
            for df in dump_files:
                dumps.append((df, LammpsDump.from_file(os.path.join(path, df))))

        # log
        log = LammpsLog(log_file=log_file)

        logger.info("Getting task doc for base dir :{}".format(path))
        d = self.generate_doc(path, lmps_input, log, dumps)

        lmps_run = None
        if len(dump_files) == 1 and data_file:
            lmps_run = LammpsRun(data_file, dump_files[0], log_file=log_filename)

        self.post_process(d, lmps_run)

        return d

    def post_process(self, d, lmps_run):
        """
        Simple post processing.

        Args:
            d (dict)
            lmps_run (LammpsRun)
        """
        if self.diffusion_params and isinstance(lmps_run, LammpsRun):
            d["analysis"] = {}
            d["analysis"]["diffusion_params"] = self.diffusion_params
            diffusion_analyzer = lmps_run.get_diffusion_analyzer(**self.diffusion_params)
            d["analysis"]["diffusion"] = diffusion_analyzer.get_summary_dict()
        d['state'] = 'successful'

    def generate_doc(self, dir_name, lmps_input, log, dumps):
        """

        Args:
            dir_name (str): path to the run dir.
            lmps_input (LammpsInput/LammpsInputSet):
            log (LammpsLog):
            dumps ([(filename, LammpsDump)]): list of (dump filename, dump) tuples

        Returns:
            dict
        """
        try:
            fullpath = os.path.abspath(dir_name)
            if self.use_full_uri:
                fullpath = get_uri(dir_name)
            d = {k: v for k, v in self.additional_fields.items()}
            d["schema"] = {"code": "atomate", "version": LammpsDrone.__version__}
            d["completed_at"] = str(datetime.fromtimestamp(os.path.getmtime(log.log_file)))
            d["dir_name"] = fullpath
            d["last_updated"] = datetime.utcnow()
            d["input"] = lmps_input.as_dict()
            d["output"] = {"log": log.as_dict()}
            d["output"]["dumps"] = dict([(dump_fname, dmp.as_dict()) for dump_fname, dmp in dumps])
            return d

        except:
            import traceback
            logger.error(traceback.format_exc())
            logger.error("Error in " + os.path.abspath(dir_name) + ".\n" + traceback.format_exc())
            return None

    def get_valid_paths(self, path):
        return [path]

    def as_dict(self):
        init_args = {"additional_fields": self.additional_fields,
                     "use_full_uri": self.use_full_uri,
                     "diffusion_params": self.diffusion_params}

        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "version": self.__class__.__version__,
                "init_args": init_args}

    @classmethod
    def from_dict(cls, d):
        return cls(**d["init_args"])
