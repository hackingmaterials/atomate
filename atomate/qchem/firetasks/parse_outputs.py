import json
import os
import shutil
import struct
from monty.io import zopen

from fireworks import FiretaskBase, FWAction, explicit_serialize
from fireworks.utilities.fw_serializers import DATETIME_HANDLER

from atomate.common.firetasks.glue_tasks import get_calc_loc
from atomate.qchem.database import QChemCalcDb
from atomate.qchem.drones import QChemDrone
from atomate.utils.utils import env_chk, get_logger

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "4/25/18"
__credits__ = "Brandon Wood, Shyam Dwaraknath, Xiaohui Qu"

logger = get_logger(__name__)


@explicit_serialize
class QChemToDb(FiretaskBase):
    """
    Enter a QChem run into the database. Uses current directory unless you
    specify calc_dir or calc_loc.

    Optional params:
        calc_dir (str): path to dir (on current filesystem) that contains QChem
            input and output files. Default: use current working directory.
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
            search for the most recent calc_loc with the matching name
        input_file (str): name of the QChem input file
        output_file (str): name of the QChem output file
        additional_fields (dict): dict of additional fields to add
        db_file (str): path to file containing the database credentials.
            Supports env_chk. Default: write data to JSON file.
        fw_spec_field (str): if set, will update the task doc with the contents
            of this key in the fw_spec.
        multirun (bool): Whether the job to parse includes multiple
            calculations in one input / output pair.
        runs (list): Series of file suffixes that the Drone should look for
            when parsing output.
    """

    optional_params = [
        "calc_dir",
        "calc_loc",
        "input_file",
        "output_file",
        "parse_grad_file",
        "parse_hess_file",
        "additional_fields",
        "db_file",
        "fw_spec_field",
        "multirun",
        "runs",
    ]

    def run_task(self, fw_spec):
        # get the directory that contains the QChem dir to parse
        calc_dir = os.getcwd()
        if "calc_dir" in self:
            calc_dir = self["calc_dir"]
        elif self.get("calc_loc"):
            calc_dir = get_calc_loc(self["calc_loc"], fw_spec["calc_locs"])["path"]
        input_file = self.get("input_file", "mol.qin")
        output_file = self.get("output_file", "mol.qout")
        multirun = self.get("multirun", False)
        runs = self.get("runs", None)

        # parse the QChem directory
        logger.info(f"PARSING DIRECTORY: {calc_dir}")

        additional_fields = self.get("additional_fields", {})

        drone = QChemDrone(runs=runs, additional_fields=additional_fields)

        # assimilate (i.e., parse)
        task_doc = drone.assimilate(
            path=calc_dir,
            input_file=input_file,
            output_file=output_file,
            multirun=multirun,
        )

        logger.info("Drone finished!")

        # parse the GRAD file, if desired and if it is present
        if self.get("parse_grad_file", False):
            grad_file = None
            if os.path.exists(os.path.join(calc_dir, "GRAD.gz")):
                grad_file = os.path.join(calc_dir, "GRAD.gz")
            elif os.path.exists(os.path.join(calc_dir, "GRAD")):
                grad_file = os.path.join(calc_dir, "GRAD")
            elif os.path.exists(os.path.join(calc_dir, "scratch/GRAD.gz")):
                grad_file = os.path.join(calc_dir, "scratch/GRAD.gz")
            elif os.path.exists(os.path.join(calc_dir, "scratch/GRAD")):
                grad_file = os.path.join(calc_dir, "scratch/GRAD")
            elif os.path.exists(os.path.join(calc_dir, "131.0.gz")):
                grad_file = os.path.join(calc_dir, "131.0.gz")
            elif os.path.exists(os.path.join(calc_dir, "131.0")):
                grad_file = os.path.join(calc_dir, "131.0")
            elif os.path.exists(os.path.join(calc_dir, "scratch/131.0.gz")):
                grad_file = os.path.join(calc_dir, "scratch/131.0.gz")
            elif os.path.exists(os.path.join(calc_dir, "scratch/131.0")):
                grad_file = os.path.join(calc_dir, "scratch/131.0")

            if grad_file is None:
                task_doc["warnings"]["grad_file_missing"] = True
            else:
                logger.info("Parsing gradient file")
                grad = []
                if grad_file[-5:] == "131.0" or grad_file[-8:] == "131.0.gz":
                    tmp_grad_data = []
                    logger.info("Gradient file is binary")
                    with zopen(grad_file, mode="rb") as file:
                        binary = file.read()
                        for ii in range(int(len(binary) / 8)):
                            tmp_grad_data.append(
                                struct.unpack("d", binary[ii * 8 : (ii + 1) * 8])[0]
                            )
                    grad = []
                    for ii in range(int(len(tmp_grad_data) / 3)):
                        grad.append(
                            [
                                float(tmp_grad_data[ii * 3]),
                                float(tmp_grad_data[ii * 3 + 1]),
                                float(tmp_grad_data[ii * 3 + 2]),
                            ]
                        )
                else:
                    logger.info("Gradient file is text")
                    with zopen(grad_file, mode="rt", encoding="ISO-8859-1") as f:
                        lines = f.readlines()
                        for line in lines:
                            split_line = line.split()
                            if len(split_line) == 3:
                                grad.append(
                                    [
                                        float(split_line[0]),
                                        float(split_line[1]),
                                        float(split_line[2]),
                                    ]
                                )
                task_doc["output"]["precise_gradients"] = grad
                logger.info("Done parsing gradient. Now removing scratch directory")
                if os.path.exists(os.path.join(calc_dir, "scratch")):
                    shutil.rmtree(os.path.join(calc_dir, "scratch"))

        # parse the HESS file(s), if desired and if one or multiple are present
        if self.get("parse_hess_file", False):
            hess_files = []
            for filename in os.listdir(calc_dir):
                if "HESS" in filename:
                    hess_files.append(filename)
                elif filename == "scratch":
                    for subfilename in os.listdir(os.path.join(calc_dir, "scratch")):
                        if "HESS" in subfilename:
                            hess_files.append("scratch/" + subfilename)
                        elif subfilename == "132.0" or subfilename == "132.0.gz":
                            hess_files.append("scratch/" + subfilename)

            if len(hess_files) == 0:
                task_doc["warnings"]["hess_file_missing"] = True
            else:
                logger.info("Parsing Hessian file")
                hess_data = {}
                for hess_name in hess_files:
                    if hess_name[-5:] == "132.0" or hess_name[-8:] == "132.0.gz":
                        logger.info("Hessian file is binary")
                        tmp_hess_data = []
                        with zopen(
                            os.path.join(calc_dir, hess_name), mode="rb"
                        ) as file:
                            binary = file.read()
                            for ii in range(int(len(binary) / 8)):
                                tmp_hess_data.append(
                                    struct.unpack("d", binary[ii * 8 : (ii + 1) * 8])[0]
                                )
                        hess_data[hess_name] = tmp_hess_data
                    else:
                        logger.info("Hessian file is text")
                        with zopen(
                            os.path.join(calc_dir, hess_name),
                            mode="rt",
                            encoding="ISO-8859-1",
                        ) as f:
                            hess_data[hess_name] = f.readlines()
                task_doc["output"]["hess_data"] = hess_data
                logger.info("Done parsing Hessian. Now removing scratch directory")
                if os.path.exists(os.path.join(calc_dir, "scratch")):
                    shutil.rmtree(os.path.join(calc_dir, "scratch"))

        # Check for additional keys to set based on the fw_spec
        logger.info("Updating fw_spec_field")
        if self.get("fw_spec_field"):
            task_doc.update(
                {self.get("fw_spec_field"): fw_spec.get(self.get("fw_spec_field"))}
            )

        # Update fw_spec with final/optimized structure
        logger.info("Updating charges in spec")
        update_spec = {}
        if task_doc.get("output").get("optimized_molecule"):
            update_spec["prev_calc_molecule"] = task_doc["output"]["optimized_molecule"]
            update_spec["prev_calc_mulliken"] = task_doc["output"]["mulliken"]
            if "RESP" in task_doc["output"]:
                update_spec["prev_calc_resp"] = task_doc["output"]["RESP"]
            elif "ESP" in task_doc["output"]:
                update_spec["prev_calc_esp"] = task_doc["output"]["ESP"]

        # get the database connection
        logger.info("Get database connection")
        db_file = env_chk(self.get("db_file"), fw_spec)

        # db insertion or taskdoc dump
        if not db_file:
            with open(os.path.join(calc_dir, "task.json"), "w") as f:
                f.write(json.dumps(task_doc, default=DATETIME_HANDLER))
        else:
            logger.info("Connecting to QChemCalcDb")
            mmdb = QChemCalcDb.from_db_file(db_file, admin=True)
            logger.info("Starting task_doc insert")
            t_id = mmdb.insert(task_doc)
            logger.info(f"Finished parsing with task_id: {t_id}")

        return FWAction(
            stored_data={"task_id": task_doc.get("task_id", None)},
            update_spec=update_spec,
        )


@explicit_serialize
class ProtCalcToDb(FiretaskBase):
    """
    Enter a QChem run of a proton into the database. Uses current directory unless you
    specify calc_dir or calc_loc.

    Optional params:
        calc_dir (str): path to dir (on current filesystem) that contains QChem
            input and output files for both the H2O and H3O+ calculations. Default: use current working directory.
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
            search for the most recent calc_loc with the matching name
        input_file_H2O (str): name of the QChem input file for H2O calculation
        output_file_H2O (str): name of the QChem output file for H2O calculation
        input_file_H3O (str): name of the QChem input file for H3O+ calculation
        output_file_H3O (str): name of the QChem output file for H3O+ calculation
        additional_fields (dict): dict of additional fields to add
        db_file (str): path to file containing the database credentials.
            Supports env_chk. Default: write data to JSON file.
        runs (list): Series of file suffixes that the Drone should look for
            when parsing output.
    """

    optional_params = [
        "calc_dir",
        "calc_loc",
        "input_file_H2O",
        "output_file_H2O",
        "input_file_H3O",
        "output_file_H3O",
        "additional_fields",
        "db_file",
        "runs",
    ]

    def run_task(self, fw_spec):
        # get the directory that contains the QChem dir to parse
        calc_dir = os.getcwd()
        if "calc_dir" in self:
            calc_dir = self["calc_dir"]
        elif self.get("calc_loc"):
            calc_dir = get_calc_loc(self["calc_loc"], fw_spec["calc_locs"])["path"]
        input_file_H2O = self.get("input_file_H2O", "water.qin")
        output_file_H2O = self.get("output_file_H2O", "water.qout")
        input_file_H3O = self.get("input_file_H3O", "hydronium.qin")
        output_file_H3O = self.get("output_file_H3O", "hydronium.qout")
        runs = self.get("runs", None)

        # parse the QChem directory
        logger.info(f"PARSING DIRECTORY: {calc_dir}")

        additional_fields = self.get("additional_fields", {})

        drone = QChemDrone(runs=runs, additional_fields=additional_fields)

        # assimilate (i.e., parse)
        task_doc_1 = drone.assimilate(
            path=calc_dir,
            input_file=input_file_H2O,
            output_file=output_file_H2O,
            multirun=False,
        )

        task_doc_2 = drone.assimilate(
            path=calc_dir,
            input_file=input_file_H3O,
            output_file=output_file_H3O,
            multirun=False,
        )

        task_doc_clean = task_doc_1
        task_doc_clean["calcs_reversed"].append(task_doc_2["calcs_reversed"][0])
        task_doc_clean["input"]["initial_molecule"]["charge"] = 1
        task_doc_clean["input"]["initial_molecule"]["spin_multiplicity"] = 1
        task_doc_clean["orig"]["molecule"]["charge"] = 1
        task_doc_clean["orig"]["molecule"]["spin_multiplicity"] = 1
        task_doc_clean["output"]["initial_molecule"]["charge"] = 1
        task_doc_clean["output"]["initial_molecule"]["spin_multiplicity"] = 1
        task_doc_clean["output"]["initial_molecule"]["sites"] = [
            {
                "name": "H",
                "species": [{"element": "H", "occu": 1}],
                "xyz": [0.0, 0.0, 0.0],
                "properties": {},
            }
        ]
        task_doc_clean["output"]["mulliken"] = [+1.000000]
        task_doc_clean["output"]["resp"] = [+1.000000]
        task_doc_clean["output"]["optimized_molecule"]["charge"] = 1
        task_doc_clean["output"]["optimized_molecule"]["spin_multiplicity"] = 1
        task_doc_clean["output"]["optimized_molecule"]["sites"] = [
            {
                "name": "H",
                "species": [{"element": "H", "occu": 1}],
                "xyz": [0.0, 0.0, 0.0],
                "properties": {},
            }
        ]
        task_doc_clean["output"]["final_energy"] = (
            task_doc_2["output"]["final_energy"] - task_doc_1["output"]["final_energy"]
        )

        # get the database connection
        db_file = env_chk(self.get("db_file"), fw_spec)

        # db insertion or taskdoc dump
        if not db_file:
            with open(os.path.join(calc_dir, "task.json"), "w") as f:
                f.write(json.dumps(task_doc_clean, default=DATETIME_HANDLER))
        else:
            mmdb = QChemCalcDb.from_db_file(db_file, admin=True)
            t_id = mmdb.insert(task_doc_clean)
            logger.info(f"Finished parsing with task_id: {t_id}")

        return FWAction(
            stored_data={"task_id": task_doc_clean.get("task_id", None)},
        )
