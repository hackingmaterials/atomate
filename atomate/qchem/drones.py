# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import datetime
from fnmatch import fnmatch
from collections import OrderedDict
import json
import glob
import traceback
from itertools import chain

from monty.io import zopen
from monty.json import jsanitize
from pymatgen.io.qchem.outputs import QCOutput
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.apps.borg.hive import AbstractDrone
from pymatgen.io.babel import BabelMolAdaptor
from pymatgen.symmetry.analyzer import PointGroupAnalyzer

from atomate.utils.utils import get_logger
from atomate import __version__ as atomate_version

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "4/25/18"
__credits__ = "Brandon Wood, Shyam Dwaraknath, Xiaohui Qu, Kiran Mathew, Shyue Ping Ong, Anubhav Jain"

logger = get_logger(__name__)


class QChemDrone(AbstractDrone):
    """
    A QChem drone to parse QChem calculations and insert an organized, searchable entry into the database.
    """

    __version__ = atomate_version  # note: the version is inserted into the task doc

    # Schema def of important keys and sub-keys; used in validation
    schema = {
        "root": {
            "dir_name", "input", "output", "calcs_reversed", "smiles",
            "walltime", "cputime", "formula_pretty", "formula_anonymous",
            "chemsys", "pointgroup"
        },
        "input": {"initial_molecule", "job_type"},
        "output": {"initial_molecule", "job_type"}
    }

    def __init__(self, runs=None, additional_fields=None):
        """
        Initialize a QChem drone to parse qchem calculations
        Args:
            runs (list): Naming scheme for multiple calcuations in one folder
            additional_fields (dict): dictionary of additional fields to add to output document
        """
        self.runs = runs or list(
            chain.from_iterable([["opt_" + str(ii), "freq_" + str(ii)]
                                 for ii in range(9)]))
        self.additional_fields = additional_fields or {}

    def assimilate(self, path, input_file, output_file, multirun):
        """
        Parses qchem input and output files and insert the result into the db.

        Args:
            path (str): Path to the directory containing output file
            input_file (str): base name of the input file(s)
            output_file (str): base name of the output file(s)
            multirun (bool): Whether the job to parse includes multiple
                             calculations in one input / output pair.

        Returns:
            d (dict): a task dictionary
        """
        logger.info("Getting task doc for base dir :{}".format(path))
        qcinput_files = self.filter_files(path, file_pattern=input_file)
        qcoutput_files = self.filter_files(path, file_pattern=output_file)
        if len(qcinput_files) != len(qcoutput_files):
            raise AssertionError("Inequal number of input and output files!")
        if len(qcinput_files) > 0 and len(qcoutput_files) > 0:
            d = self.generate_doc(path, qcinput_files, qcoutput_files,
                                  multirun)
            self.post_process(path, d)
        else:
            raise ValueError("Either input or output not found!")
        self.validate_doc(d)
        return jsanitize(d, strict=True, allow_bson=True)

    def filter_files(self, path, file_pattern):
        """
        Find the files that match the pattern in the given path and
        return them in an ordered dictionary. The searched for files are
        filtered by the run types defined in self.runs.

        Args:
            path (string): path to the folder
            file_pattern (string): base files to be searched for

        Returns:
            OrderedDict of the names of the files to be processed further.
            The key is set from list of run types: self.runs
        """
        processed_files = OrderedDict()
        files = os.listdir(path)
        for r in self.runs:
            # try subfolder schema
            if r in files:
                for f in os.listdir(os.path.join(path, r)):
                    if fnmatch(f, "{}*".format(file_pattern)):
                        processed_files[r] = os.path.join(r, f)
            # try extension schema
            else:
                for f in files:
                    if fnmatch(f, "{}.{}*".format(file_pattern, r)):
                        processed_files[r] = f
        if len(processed_files) == 0:
            # get any matching file from the folder
            for f in files:
                if fnmatch(f, "{}*".format(file_pattern)):
                    processed_files["standard"] = f
        return processed_files

    def generate_doc(self, dir_name, qcinput_files, qcoutput_files, multirun):
        try:
            fullpath = os.path.abspath(dir_name)
            d = jsanitize(self.additional_fields, strict=True)
            d["schema"] = {
                "code": "atomate",
                "version": QChemDrone.__version__
            }
            d["dir_name"] = fullpath
            if multirun:
                d["calcs_reversed"] = self.process_qchem_multirun(
                    dir_name, qcinput_files, qcoutput_files)
            else:
                d["calcs_reversed"] = [
                    self.process_qchemrun(dir_name, taskname,
                                          qcinput_files.get(taskname),
                                          output_filename)
                    for taskname, output_filename in qcoutput_files.items()
                ]

            # reverse the calculations data order so newest calc is first
            d["calcs_reversed"].reverse()

            d_calc_init = d["calcs_reversed"][-1]
            d_calc_final = d["calcs_reversed"][0]

            d["input"] = {
                "initial_molecule": d_calc_init["initial_molecule"],
                "job_type": d_calc_init["input"]["rem"]["job_type"]
            }
            d["output"] = {
                "initial_molecule": d_calc_final["initial_molecule"],
                "job_type": d_calc_final["input"]["rem"]["job_type"]
            }

            if d["output"]["job_type"] == "opt" or d["output"]["job_type"] == "optimization":
                d["output"]["optimized_molecule"] = d_calc_final[
                    "molecule_from_optimized_geometry"]
                d["output"]["final_energy"] = d_calc_final["final_energy"]
                if d_calc_final["opt_constraint"]:
                    d["output"]["constraint"] = [
                        d_calc_final["opt_constraint"][0],
                        float(d_calc_final["opt_constraint"][6])
                    ]
            if d["output"]["job_type"] == "freq" or d["output"]["job_type"] == "frequency":
                d["output"]["frequencies"] = d_calc_final["frequencies"]
                d["output"]["enthalpy"] = d_calc_final["enthalpy"]
                d["output"]["entropy"] = d_calc_final["entropy"]
                if d["input"]["job_type"] == "opt" or d["input"]["job_type"] == "optimization":
                    d["output"]["optimized_molecule"] = d_calc_final[
                        "initial_molecule"]
                    d["output"]["final_energy"] = d["calcs_reversed"][1][
                        "final_energy"]

            if "special_run_type" in d:
                if d["special_run_type"] == "frequency_flattener":
                    d["num_frequencies_flattened"] = (
                        len(qcinput_files) / 2) - 1

            total_cputime = 0.0
            total_walltime = 0.0
            nan_found = False
            for calc in d["calcs_reversed"]:
                if calc["walltime"] != "nan":
                    total_walltime += calc["walltime"]
                else:
                    nan_found = True
                if calc["cputime"] != "nan":
                    total_cputime += calc["cputime"]
                else:
                    nan_found = True
            if nan_found:
                d["walltime"] = "nan"
                d["cputime"] = "nan"
            else:
                d["walltime"] = total_walltime
                d["cputime"] = total_cputime

            comp = d["output"]["initial_molecule"].composition
            d["formula_pretty"] = comp.reduced_formula
            d["formula_anonymous"] = comp.anonymized_formula
            d["chemsys"] = "-".join(sorted(set(d_calc_final["species"])))
            d["pointgroup"] = PointGroupAnalyzer(
                d["output"]["initial_molecule"]).sch_symbol

            bb = BabelMolAdaptor(d["output"]["initial_molecule"])
            pbmol = bb.pybel_mol
            smiles = pbmol.write(str("smi")).split()[0]
            d["smiles"] = smiles

            d["state"] = "successful" if d_calc_final[
                "completion"] else "unsuccessful"
            d["last_updated"] = datetime.datetime.utcnow()
            return d

        except Exception:
            logger.error(traceback.format_exc())
            logger.error("Error in " + os.path.abspath(dir_name) + ".\n" +
                         traceback.format_exc())
            raise

    @staticmethod
    def process_qchemrun(dir_name, taskname, input_file, output_file):
        """
        Process a QChem calculation, aka an input/output pair.
        """
        qchem_input_file = os.path.join(dir_name, input_file)
        qchem_output_file = os.path.join(dir_name, output_file)
        d = QCOutput(qchem_output_file).data
        temp_input = QCInput.from_file(qchem_input_file)
        d["input"] = {}
        d["input"]["molecule"] = temp_input.molecule
        d["input"]["rem"] = temp_input.rem
        d["input"]["opt"] = temp_input.opt
        d["input"]["pcm"] = temp_input.pcm
        d["input"]["solvent"] = temp_input.solvent
        d["task"] = {"type": taskname, "name": taskname}
        return d

    @staticmethod
    def process_qchem_multirun(dir_name, input_files, output_files):
        """
        Process a QChem run which is known to include multiple calculations
        in a single input/output pair.
        """
        if len(input_files) != 1:
            raise ValueError(
                "ERROR: The drone can only process a directory containing a single input/output pair when each include multiple calculations."
            )
        else:
            for key in input_files:
                to_return = []
                qchem_input_file = os.path.join(dir_name, input_files.get(key))
                qchem_output_file = os.path.join(dir_name,
                                                 output_files.get(key))
                multi_out = QCOutput.multiple_outputs_from_file(
                    QCOutput, qchem_output_file, keep_sub_files=False)
                multi_in = QCInput.from_multi_jobs_file(qchem_input_file)
                for ii, out in enumerate(multi_out):
                    d = out.data
                    d["input"] = {}
                    d["input"]["molecule"] = multi_in[ii].molecule
                    d["input"]["rem"] = multi_in[ii].rem
                    d["input"]["opt"] = multi_in[ii].opt
                    d["input"]["pcm"] = multi_in[ii].pcm
                    d["input"]["solvent"] = multi_in[ii].solvent
                    d["task"] = {"type": key, "name": "calc" + str(ii)}
                    to_return.append(d)
            return to_return

    @staticmethod
    def post_process(dir_name, d):
        """
        Post-processing for various files other than the QChem input and output files.
        Currently only looks for custodian.json. Modify this if other files need to be processed.
        """
        logger.info("Post-processing dir:{}".format(dir_name))
        fullpath = os.path.abspath(dir_name)
        filenames = glob.glob(os.path.join(fullpath, "custodian.json*"))
        if len(filenames) >= 1:
            with zopen(filenames[0], "rt") as f:
                d["custodian"] = json.load(f)

    def validate_doc(self, d):
        """
        Sanity check, aka make sure all the important keys are set. Note that a failure
        to pass validation is unfortunately unlikely to be noticed by a user.
        """
        for k, v in self.schema.items():
            diff = v.difference(set(d.get(k, d).keys()))
            if diff:
                logger.warn("The keys {0} in {1} not set".format(diff, k))

    @staticmethod
    def get_valid_paths(self, path):
        return [path]
