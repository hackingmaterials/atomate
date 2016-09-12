# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import json
import os

from pymatgen.io.feff.outputs import Xmu

from fireworks import FireTaskBase, FWAction, explicit_serialize
from fireworks.utilities.fw_serializers import DATETIME_HANDLER

from matmethods.utils.utils import env_chk, get_calc_loc
from matmethods.utils.utils import get_logger
from matmethods.vasp.database import MMDb

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


@explicit_serialize
class XmuToDbTask(FireTaskBase):
    """
    Parse the output of EXAFS calculation(xmu.dat) and insert it into the database("xas" collection).

    Required_params:
        absorbing_atom (str): absorbing atom symbol
        structure (Structure): input structure

    Optional_params:
        input_file (str): path to the feff input file.
        calc_dir (str): path to dir (on current filesystem) that contains FEFF output files.
            Default: use current working directory.
        calc_loc (str OR bool): if True will set most recent calc_loc. If str search for the most
            recent calc_loc with the matching name
        db_file (str): path to the db file.
    """

    required_params = ["absorbing_atom", "structure"]
    optional_params = ["input_file", "calc_dir", "calc_loc", "db_file"]

    def run_task(self, fw_spec):
        calc_dir = os.getcwd()
        if "calc_dir" in self:
            calc_dir = self["calc_dir"]
        elif self.get("calc_loc"):
            calc_dir = get_calc_loc(self["calc_loc"], fw_spec["calc_locs"])["path"]

        logger.info("PARSING DIRECTORY: {}".format(calc_dir))

        db_file = env_chk(self.get('db_file'), fw_spec)

        doc = {"structure": self["structure"].as_dict(), "absorbing_atom": self["absorbing_atom"]}
        input_file = os.path.join(calc_dir, self.get("input_file", "feff.inp"))
        xmu = Xmu.from_file(os.path.join(calc_dir, "xmu.dat"), input_file)
        doc["xmu"] = xmu.data.tolist()

        # db insertion
        if not db_file:
            with open("xmu.json", "w") as f:
                f.write(json.dumps(doc, default=DATETIME_HANDLER))
        else:
            db = MMDb.from_db_file(db_file, admin=True)
            db.collection = db.db["xas"]
            db.collection.insert_one(doc)

        logger.info("Finished parsing EXAFS spectrum")

        return FWAction()
