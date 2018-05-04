# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import re
import datetime
from fnmatch import fnmatch
from collections import OrderedDict
import json # using!
import glob # using!
import traceback

from monty.io import zopen
from monty.json import jsanitize

import numpy as np

from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.core.operations import SymmOp
from pymatgen.io.qchem_io.outputs import QCOutput
from pymatgen.io.qchem_io.inputs import QCInput
from pymatgen.apps.borg.hive import AbstractDrone
from pymatgen.io.babel import BabelMolAdaptor

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
        "root": {"dir_name", "input", "output", "calcs_reversed", "smiles"},
        "input": {"initial_molecule", "job_type"},
        "output": {"initial_molecule", "job_type"}
    }

    def __init__(self, runs, additional_fields=None):
        """
        Initialize a QChem drone to parse qchem calculations
        Args:
            runs (list): Naming scheme for multiple calcuations in one folder
            additional_fields (dict): dictionary of additional fields to add to output document
        """
        self.runs = runs
        self.additional_fields = additional_fields or {}


    def assimilate(self, path, input_file, output_file):
        """
        Parses qchem input and output files and insert the result into the db.

        Args:
            path (str): Path to the directory containing output file
            input_file (str): base name of the input file(s)
            output_file (str): base name of the output file(s)

        Returns:
            d (dict): a task dictionary
        """
        logger.info("Getting task doc for base dir :{}".format(path))
        qcinput_files = self.filter_files(path, file_pattern=input_file)
        qcoutput_files = self.filter_files(path, file_pattern=output_file)
        if len(qcinput_files) != len(qcoutput_files):
            raise AssertionError('Inequal number of input and output files!')
        if len(qcinput_files) > 0 and len(qcoutput_files) > 0:
            d = self.generate_doc(path, qcinput_files, qcoutput_files)
            self.post_process(path, d)
        else:
            raise ValueError("Either input or output not found!")
        self.validate_doc(d)
        return d


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
                    processed_files['standard'] = f
        return processed_files


    def generate_doc(self, dir_name, qcinput_files, qcoutput_files):
        try:

            fullpath = os.path.abspath(dir_name)
            d = jsanitize(self.additional_fields, strict=True)
            d["schema"] = {"code": "atomate", "version": QChemDrone.__version__}
            d["dir_name"] = fullpath
            d["calcs_reversed"] = [self.process_qchemrun(dir_name, taskname, qcinput_files.get(taskname), output_filename)
                                   for taskname, output_filename in qcoutput_files.items()]

            # reverse the calculations data order so newest calc is first
            d["calcs_reversed"].reverse()

            d_calc_init = d["calcs_reversed"][-1]
            d_calc_final = d["calcs_reversed"][0]

            d["input"] = {"initial_molecule": d_calc_init["initial_molecule"],
                          "job_type": d_calc_init["input_rem"]["job_type"]}
            d["output"] = {"initial_molecule": d_calc_final["initial_molecule"],
                          "job_type": d_calc_final["input_rem"]["job_type"]}

            if d["output"]["job_type"] == "opt" or d["output"]["job_type"] == "optimization":
                d["output"]["optimized_molecule"] = d_calc_final["molecule_from_optimized_geometry"]
            if d["output"]["job_type"] == "freq" or d["output"]["job_type"] == "frequency":
                d["output"]["frequencies"] = d_calc_final["frequencies"]
                if d["input"]["job_type"] == "opt" or d["input"]["job_type"] == "optimization":
                    d["output"]["optimized_molecule"] = d_calc_final["initial_molecule"]

            if "special_run_type" in d:
                if d["special_run_type"] == "frequency_flattener":
                    d["num_frequencies_flattened"] = (len(qcinput_files)/2)-1

            bb = BabelMolAdaptor(d["output"]["initial_molecule"])
            pbmol = bb.pybel_mol
            smiles = pbmol.write(str("smi")).split()[0]
            d["smiles"] = smiles

            d["last_updated"] = datetime.datetime.utcnow()
            return d

        except Exception:
            logger.error(traceback.format_exc())
            logger.error("Error in " + os.path.abspath(dir_name) + ".\n" + traceback.format_exc())
            raise

    def process_qchemrun(self, dir_name, taskname, input_file, output_file):
        """
        Process a QChem calculation, aka an input/output pair.
        """
        qchem_input_file = os.path.join(dir_name, input_file)
        qchem_output_file = os.path.join(dir_name, output_file)
        d = QCOutput(qchem_output_file).data
        temp_input = QCInput.from_file(qchem_input_file)
        d["input_molecule"] = temp_input.molecule
        d["input_rem"] = temp_input.rem
        d["input_opt"] = temp_input.opt
        d["input_pcm"] = temp_input.pcm
        d["input_solvent"] = temp_input.solvent
        d["task"] = {"type": taskname, "name": taskname}
        return d

    def post_process(self, dir_name, d):
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

    def get_valid_paths(self, path):
        return [path]
