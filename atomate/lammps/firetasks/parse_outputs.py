# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines firetasks for parsing and processing the LAMMPS outputfiles to extract useful
information such as the summary of transport properties and insert them into the database.
"""

import json
import os

from datetime import datetime

from fireworks import FiretaskBase, FWAction
from fireworks.utilities.fw_serializers import DATETIME_HANDLER
from fireworks.utilities.fw_utilities import explicit_serialize

from atomate.utils.utils import get_logger
from atomate.common.firetasks.glue_tasks import get_calc_loc
from atomate.utils.utils import env_chk
from atomate.lammps.database import LammpsCalcDb
from pymatgen.io.lammps.output import LammpsRun


__author__ = 'Kiran Mathew'
__email__ = "kmathew@lbl.gov"

logger = get_logger(__name__)

# TODO: @matk86 - for your consideration. You can implement as a Drone and simply use ToDbTask (in
# atomate.common) to do the processing instead of a special LammpsToDBTask. Advantage is that the
# parsing code is then in pymatgen where more people could use and modify it. -computron


@explicit_serialize
class LammpsToDBTask(FiretaskBase):
    """
    Enter a LAMMPS run into the database.

    Require params:
        lammps_input (DictLammpsInput)

    Optional params:
        calc_dir (str): path to dir (on current filesystem) that contains LAMMPS
            output files. Default: use current working directory.
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
            search for the most recent calc_loc with the matching name
        diffusion_params (dict): parameters to the diffusion_analyzer. If specified a summary
            of diffusion statistics will be added.
        db_file (str): path to file containing the database credentials.
            Supports env_chk. Default: write data to JSON file.
    """

    required_params = ["lammps_input"]
    optional_params = ["calc_dir", "calc_loc", "diffusion_params", "db_file"]

    def run_task(self, fw_spec):
        lammps_input = self["lammps_input"]
        diffusion_params = self.get("diffusion_params", {})

        # get the directory that contains the LAMMPS dir to parse
        calc_dir = os.getcwd()
        if "calc_dir" in self:
            calc_dir = self["calc_dir"]
        elif self.get("calc_loc"):
            calc_dir = get_calc_loc(self["calc_loc"], fw_spec["calc_locs"])["path"]

        # parse the directory
        logger.info("PARSING DIRECTORY: {}".format(calc_dir))
        d = {}
        d["dir_name"] = os.path.abspath(os.getcwd())
        d["last_updated"] = datetime.today()
        d["input"] = lammps_input.as_dict()
        log_file = lammps_input.config_dict["log"]
        if isinstance(lammps_input.config_dict["dump"], list):
            dump_file = lammps_input.config_dict["dump"][0].split()[4]
        else:
            dump_file = lammps_input.config_dict["dump"].split()[4]
        is_forcefield = hasattr(lammps_input.lammps_data, "bonds_data")
        lammpsrun = LammpsRun(lammps_input.data_filename, dump_file, log_file,
                              is_forcefield=is_forcefield)
        d["natoms"] = lammpsrun.natoms
        d["nmols"] = lammpsrun.nmols
        d["box_lengths"] = lammpsrun.box_lengths
        d["mol_masses"] = lammpsrun.mol_masses
        d["mol_config"] = lammpsrun.mol_config
        if diffusion_params:
            diffusion_analyzer = lammpsrun.get_diffusion_analyzer(**diffusion_params)
            d["analysis"]["diffusion"] = diffusion_analyzer.get_summary_dict()
        db_file = env_chk(self.get('db_file'), fw_spec)

        # db insertion
        if not db_file:
            with open("task.json", "w") as f:
                f.write(json.dumps(d, default=DATETIME_HANDLER))
        else:
            mmdb = LammpsCalcDb.from_db_file(db_file)
            # insert the task document
            t_id = mmdb.insert(d)
            logger.info("Finished parsing with task_id: {}".format(t_id))
        return FWAction(stored_data={"task_id": d.get("task_id", None)})
