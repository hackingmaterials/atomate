# coding: utf-8


import json
import os
from datetime import datetime
from glob import glob

import numpy as np

from pymatgen.io.feff.inputs import Tags, Atoms

from fireworks import FiretaskBase, FWAction, explicit_serialize
from fireworks.utilities.fw_serializers import DATETIME_HANDLER
from fireworks.user_objects.firetasks.filepad_tasks import get_fpad

from atomate.utils.utils import env_chk
from atomate.common.firetasks.glue_tasks import get_calc_loc
from atomate.utils.utils import get_logger
from atomate.feff.database import FeffCalcDb

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


@explicit_serialize
class SpectrumToDbTask(FiretaskBase):
    """
    Parse the output of absorption/core-loss spectrum calculations(xmu.dat, eels.dat) and insert it
    into the database.

    Required_params:
        absorbing_atom (str): absorbing atom symbol
        structure (Structure): input structure
        spectrum_type (str): XANES, EXAFS, ELNES, EXELFS
        output_file (str): the output file name. xmu.dat or eels.dat

    Optional_params:
        input_file (str): path to the feff input file.
        calc_dir (str): path to dir (on current filesystem) that contains FEFF output files.
            Default: use current working directory.
        calc_loc (str OR bool): if True will set most recent calc_loc. If str search for the most
            recent calc_loc with the matching name
        db_file (str): path to the db file.
        edge (str): absorption edge
        metadata (dict): meta data
    """

    required_params = ["absorbing_atom", "structure", "spectrum_type", "output_file"]
    optional_params = ["input_file", "calc_dir", "calc_loc", "db_file", "edge", "metadata"]

    def run_task(self, fw_spec):
        calc_dir = os.getcwd()
        if "calc_dir" in self:
            calc_dir = self["calc_dir"]
        elif self.get("calc_loc"):
            calc_dir = get_calc_loc(self["calc_loc"], fw_spec["calc_locs"])["path"]

        logger.info("PARSING DIRECTORY: {}".format(calc_dir))

        db_file = env_chk(self.get('db_file'), fw_spec)

        cluster_dict = None
        tags = Tags.from_file(filename="feff.inp")
        if "RECIPROCAL" not in tags:
            cluster_dict = Atoms.cluster_from_file("feff.inp").as_dict()
        doc = {"input_parameters": tags.as_dict(),
               "cluster": cluster_dict,
               "structure": self["structure"].as_dict(),
               "absorbing_atom": self["absorbing_atom"],
               "spectrum_type": self["spectrum_type"],
               "spectrum": np.loadtxt(os.path.join(calc_dir, self["output_file"])).tolist(),
               "edge": self.get("edge", None),
               "metadata": self.get("metadata", None),
               "dir_name": os.path.abspath(os.getcwd()),
               "last_updated": datetime.utcnow()}

        if not db_file:
            with open("feff_task.json", "w") as f:
                f.write(json.dumps(doc, default=DATETIME_HANDLER))

        else:
            db = FeffCalcDb.from_db_file(db_file, admin=True)
            db.insert(doc)

        logger.info("Finished parsing the spectrum")

        return FWAction(stored_data={"task_id": doc.get("task_id", None)})


@explicit_serialize
class AddPathsToFilepadTask(FiretaskBase):
    """
    Insert the scattering amplitude outputs(all feffNNNN.dat files) to gridfs using filepad.

    Optional_params:
        labels (list): list of labels to tag the inserted files. Useful for querying later.
        filepad_file (str): path to the filepad connection settings file.
        compress (bool): wether or not to compress the file contents before insertion.
        metadata (dict): metadata.
    """

    optional_params = ["labels", "filepad_file", "compress", "metadata"]

    def run_task(self, fw_spec):
        paths = glob("feff????.dat")
        fpad = get_fpad(self.get("filepad_file", None))
        labels = self.get("labels", None)
        for i, p in enumerate(paths):
            l = labels[i] if labels is not None else None
            fpad.add_file(p, label=l, metadata=self.get("metadata", None),
                          compress=self.get("compress", True))
