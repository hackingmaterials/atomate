# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import json
import os
import copy
import datetime

from fireworks import FiretaskBase, FWAction, explicit_serialize
from fireworks.utilities.fw_serializers import DATETIME_HANDLER
from pymatgen.io.babel import BabelMolAdaptor
from monty.serialization import dumpfn
from monty.json import jsanitize

from atomate.qchem.database import QChemCalcDb
from atomate.utils.utils import env_chk
from atomate.utils.utils import get_logger

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "11/13/18"

logger = get_logger(__name__)


@explicit_serialize
class GatherGeometries(FiretaskBase):
    """
    """
    required_params = ["prefix", "suffix"]
    optional_params = ["db_file"]

    def run_task(self, fw_spec):

        data = []
        for key in fw_spec:
            if len(key) > len(self["prefix"])+len(self["suffix"]):
                if key[0:len(self["prefix"])] == self["prefix"]:
                    to_add = copy.deepcopy(fw_spec[key])
                    if len(self["suffix"]) == 0:
                        to_add["index"] = key[len(self["prefix"]):]
                    else:
                        to_add["index"] = key[len(self["prefix"]):-len(self["suffix"])]
                    data.append(to_add)
        if len(data) == 0:
            raise KeyError("ERROR: fw_spec does not contain any keys with the given prefix / suffix! Exiting...")
        else:
            sorted_data = sorted(data, key=lambda k: k["energy"])

        for ii in range(3):
            print(sorted_data[ii]["energy"])

        task_doc = {}
        task_doc["sorted_data"] = sorted_data
        task_doc["dir_name"] = os.getcwd()
        task_doc["task_label"] = "gather_geoms" + self["suffix"]
        comp = sorted_data[0]["molecule"].composition
        task_doc["formula_pretty"] = comp.reduced_formula
        task_doc["formula_anonymous"] = comp.anonymized_formula
        bb = BabelMolAdaptor(sorted_data[0]["molecule"])
        pbmol = bb.pybel_mol
        smiles = pbmol.write(str("smi")).split()[0]
        task_doc["smiles"] = smiles
        task_doc["last_updated"] = datetime.datetime.utcnow()
        task_doc = jsanitize(task_doc, strict=True, allow_bson=True)

        db_file = env_chk(self.get("db_file"), fw_spec)

        if not db_file:
            raise AssertionError("ERROR: No db file found! Dumped task doc to json. Exiting...")
        else:
            mmdb = QChemCalcDb.from_db_file(db_file, admin=True)
            t_id = mmdb.insert(task_doc)
            logger.info("Finished parsing with task_id: {}".format(t_id))



