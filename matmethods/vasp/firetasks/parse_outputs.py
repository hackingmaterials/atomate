# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import json
import os

from datetime import datetime

from monty.json import MontyEncoder

from fireworks import FireTaskBase, FWAction
from fireworks.utilities.fw_serializers import DATETIME_HANDLER
from fireworks.utilities.fw_utilities import explicit_serialize

from matmethods.utils.utils import env_chk, get_calc_loc, get_meta_from_structure
from matmethods.utils.utils import get_logger
from matmethods.vasp.database import MMDb
from matmethods.vasp.drones import VaspDrone

from pymatgen import Structure
from pymatgen.electronic_structure.boltztrap import BoltztrapAnalyzer
from pymatgen.io.vasp.sets import get_vasprun_outcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

__author__ = 'Anubhav Jain, Kiran Mathew, Shyam Dwaraknath'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov, shyamd@lbl.gov'

logger = get_logger(__name__)


@explicit_serialize
class VaspToDbTask(FireTaskBase):
    """
    Enter a VASP run into the database. Uses current directory unless you
    specify calc_dir or calc_loc.

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
        fw_spec_field (str): if set, will update the task doc with the contents
            of this key in the fw_spec.
        db_file (str): path to file containing the database credentials.
            Supports env_chk. Default: write data to JSON file.
        defuse_unsuccessful (bool): Defuses children fireworks if VASP run state
            is not "successful"; i.e. both electronic and ionic convergence are reached.
            Defaults to True.
    """
    optional_params = ["calc_dir", "calc_loc", "parse_dos", "bandstructure_mode",
                       "additional_fields", "db_file", "fw_spec_fields", "defuse_unsuccessful"]

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

        # assimilate (i.e., parse)
        task_doc = drone.assimilate(calc_dir)

        # Check for additional fields to add in the fw_spec
        if self.get("fw_spec_field"):
            task_doc.update(fw_spec[self.get("fw_spec_field")])

        # db insertion
        if not db_file:
            with open("task.json", "w") as f:
                f.write(json.dumps(task_doc, default=DATETIME_HANDLER))
        else:
            mmdb = MMDb.from_db_file(db_file, admin=True)

            # insert dos into GridFS
            if self.get("parse_dos") and "calcs_reversed" in task_doc:
                for idx, x in enumerate(task_doc["calcs_reversed"]):
                    if "dos" in task_doc["calcs_reversed"][idx]:
                        if idx == 0:  # only store most recent DOS
                            dos = json.dumps(task_doc["calcs_reversed"][idx]["dos"], cls=MontyEncoder)
                            gfs_id, compression_type = mmdb.insert_gridfs(dos, "dos_fs")
                            task_doc["calcs_reversed"][idx]["dos_compression"] = compression_type
                            task_doc["calcs_reversed"][idx]["dos_fs_id"] = gfs_id
                        del task_doc["calcs_reversed"][idx]["dos"]

            # insert band structure into GridFS
            if self.get("bandstructure_mode") and "calcs_reversed" in task_doc:
                for idx, x in enumerate(task_doc["calcs_reversed"]):
                    if "bandstructure" in task_doc["calcs_reversed"][idx]:
                        if idx == 0:  # only store most recent band structure
                            bs = json.dumps(task_doc["calcs_reversed"][idx]["bandstructure"], cls=MontyEncoder)
                            gfs_id, compression_type = mmdb.insert_gridfs(bs, "bandstructure_fs")
                            task_doc["calcs_reversed"][idx]["bandstructure_compression"] = compression_type
                            task_doc["calcs_reversed"][idx]["bandstructure_fs_id"] = gfs_id
                        del task_doc["calcs_reversed"][idx]["bandstructure"]

            # insert the task document
            t_id = mmdb.insert(task_doc)

            logger.info("Finished parsing with task_id: {}".format(t_id))

        if self.get("defuse_unsuccessful", True):
            defuse_children = (task_doc["state"] != "successful")
        else:
            defuse_children = False

        return FWAction(stored_data={"task_id": task_doc.get("task_id", None)},
                        defuse_children=defuse_children)


@explicit_serialize
class BoltztrapToDBTask(FireTaskBase):
    """
    Enter a Boltztrap run into the database.

    Optional params:
        db_file (str): path to file containing the database credentials.
            Supports env_chk. Default: write data to JSON file.
        hall_doping (bool): set True to retain hall_doping in dict
    """

    optional_params = ["db_file", "hall_doping"]

    def run_task(self, fw_spec):
        btrap_dir = os.path.join(os.getcwd(), "boltztrap")
        bta = BoltztrapAnalyzer.from_files(btrap_dir)
        d = bta.as_dict()
        d["boltztrap_dir"] = btrap_dir

        # trim the output
        for x in ['cond', 'seebeck', 'kappa', 'hall', 'mu_steps',
                  'mu_doping', 'carrier_conc']:
            del d[x]

        if not self.get("hall_doping"):
            del d["hall_doping"]

        # add the structure
        bandstructure_dir = os.getcwd()
        v, o = get_vasprun_outcar(bandstructure_dir, parse_eigen=False,
                                  parse_dos=False)
        structure = v.final_structure
        d["structure"] = structure.as_dict()
        d.update(get_meta_from_structure(structure))
        d["bandstructure_dir"] = bandstructure_dir

        # add the spacegroup
        sg = SpacegroupAnalyzer(Structure.from_dict(d["structure"]), 0.1)
        d["spacegroup"] = {"symbol": sg.get_spacegroup_symbol(),
                           "number": sg.get_spacegroup_number(),
                           "point_group": sg.get_point_group(),
                           "source": "spglib",
                           "crystal_system": sg.get_crystal_system(),
                           "hall": sg.get_hall()}

        d["created_at"] = datetime.utcnow()

        db_file = env_chk(self.get('db_file'), fw_spec)

        if not db_file:
            with open(os.path.join(btrap_dir, "boltztrap.json"), "w") as f:
                f.write(json.dumps(d, default=DATETIME_HANDLER))
        else:
            mmdb = MMDb.from_db_file(db_file, admin=True)

            # dos gets inserted into GridFS
            dos = json.dumps(d["dos"], cls=MontyEncoder)
            fsid, compression = mmdb.insert_gridfs(
                dos, collection="dos_boltztrap_fs", compress=True)
            d["dos_boltztrap_fs_id"] = fsid
            del d["dos"]

            mmdb.db.boltztrap.insert(d)
