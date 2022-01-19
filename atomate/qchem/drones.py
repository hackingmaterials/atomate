import copy
import datetime
import glob
import json
import os
import traceback
from collections import OrderedDict
from fnmatch import fnmatch
from itertools import chain

from monty.io import zopen
from monty.json import jsanitize
from pymatgen.apps.borg.hive import AbstractDrone
from pymatgen.core import Molecule
from pymatgen.io.babel import BabelMolAdaptor
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.outputs import QCOutput
from pymatgen.symmetry.analyzer import PointGroupAnalyzer

from atomate import __version__ as atomate_version
from atomate.utils.utils import get_logger

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
            "dir_name",
            "input",
            "output",
            "calcs_reversed",
            "smiles",
            "walltime",
            "cputime",
            "formula_pretty",
            "formula_anonymous",
            "chemsys",
            "pointgroup",
            "formula_alphabetical",
        },
        "input": {"initial_molecule", "job_type"},
        "output": {"initial_molecule", "job_type", "final_energy"},
    }

    def __init__(self, runs=None, additional_fields=None):
        """
        Initialize a QChem drone to parse qchem calculations
        Args:
            runs (list): Naming scheme for multiple calculations in one folder
            additional_fields (dict): dictionary of additional fields to add to output document
        """
        self.runs = runs or list(
            chain.from_iterable(
                [["opt_" + str(ii), "freq_" + str(ii)] for ii in range(10)]
            )
        )
        self.runs = ["orig"] + self.runs
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
        logger.info(f"Getting task doc for base dir :{path}")
        qcinput_files = self.filter_files(path, file_pattern=input_file)
        qcoutput_files = self.filter_files(path, file_pattern=output_file)
        if len(qcinput_files) != len(qcoutput_files):
            if len(qcinput_files) > len(qcoutput_files):
                if list(qcinput_files.items())[0][0] != "orig":
                    raise AssertionError(
                        "Can only have unequal number of input and output files when there is a saved copy "
                        "of the original input!"
                    )
            else:
                raise AssertionError("Inequal number of input and output files!")
        if len(qcinput_files) > 0 and len(qcoutput_files) > 0:
            d = self.generate_doc(path, qcinput_files, qcoutput_files, multirun)
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
                    if fnmatch(f, f"{file_pattern}*"):
                        processed_files[r] = os.path.join(r, f)
            # try extension schema
            else:
                for f in files:
                    if fnmatch(f, f"{file_pattern}.{r}*"):
                        processed_files[r] = f
        if len(processed_files) == 0 or (
            len(processed_files) == 1 and "orig" in processed_files
        ):
            # get any matching file from the folder
            for f in files:
                if fnmatch(f, f"{file_pattern}*"):
                    processed_files["standard"] = f
        return processed_files

    def generate_doc(self, dir_name, qcinput_files, qcoutput_files, multirun):
        try:
            fullpath = os.path.abspath(dir_name)
            d = jsanitize(self.additional_fields, strict=True)
            d["schema"] = {"code": "atomate", "version": QChemDrone.__version__}
            d["dir_name"] = fullpath

            # If a saved "orig" input file is present, parse it in case the error handler made changes
            # to the initial input molecule or rem params, which we might want to filter for later
            if len(qcinput_files) > len(qcoutput_files):
                orig_input = QCInput.from_file(
                    os.path.join(dir_name, qcinput_files.pop("orig"))
                )
                d["orig"] = {}
                d["orig"]["molecule"] = orig_input.molecule.as_dict()
                d["orig"]["molecule"]["charge"] = int(d["orig"]["molecule"]["charge"])
                d["orig"]["rem"] = orig_input.rem
                d["orig"]["opt"] = orig_input.opt
                d["orig"]["pcm"] = orig_input.pcm
                d["orig"]["solvent"] = orig_input.solvent
                d["orig"]["smx"] = orig_input.smx
                d["orig"]["vdw_mode"] = orig_input.vdw_mode
                d["orig"]["van_der_waals"] = orig_input.van_der_waals

            if multirun:
                d["calcs_reversed"] = self.process_qchem_multirun(
                    dir_name, qcinput_files, qcoutput_files
                )
            else:
                d["calcs_reversed"] = [
                    self.process_qchemrun(
                        dir_name, taskname, qcinput_files.get(taskname), output_filename
                    )
                    for taskname, output_filename in qcoutput_files.items()
                ]

            # reverse the calculations data order so newest calc is first
            d["calcs_reversed"].reverse()

            d["structure_change"] = []
            d["warnings"] = {}
            for entry in d["calcs_reversed"]:
                if (
                    "structure_change" in entry
                    and "structure_change" not in d["warnings"]
                ):
                    if entry["structure_change"] != "no_change":
                        d["warnings"]["structure_change"] = True
                if "structure_change" in entry:
                    d["structure_change"].append(entry["structure_change"])
                for key in entry["warnings"]:
                    if key not in d["warnings"]:
                        d["warnings"][key] = True

            d_calc_init = d["calcs_reversed"][-1]
            d_calc_final = d["calcs_reversed"][0]

            d["input"] = {
                "initial_molecule": d_calc_init["initial_molecule"],
                "job_type": d_calc_init["input"]["rem"]["job_type"],
            }
            d["output"] = {
                "initial_molecule": d_calc_final["initial_molecule"],
                "job_type": d_calc_final["input"]["rem"]["job_type"],
                "mulliken": d_calc_final["Mulliken"][-1],
            }
            if "RESP" in d_calc_final:
                d["output"]["resp"] = d_calc_final["RESP"][-1]
            elif "ESP" in d_calc_final:
                d["output"]["esp"] = d_calc_final["ESP"][-1]

            if "nbo_data" in d_calc_final:
                d["output"]["nbo"] = d_calc_final["nbo_data"]

            if d["output"]["job_type"] in ["opt", "optimization", "ts"]:
                if "molecule_from_optimized_geometry" in d_calc_final:
                    d["output"]["optimized_molecule"] = d_calc_final[
                        "molecule_from_optimized_geometry"
                    ]
                    d["output"]["final_energy"] = d_calc_final["final_energy"]
                else:
                    d["output"]["final_energy"] = "unstable"
                if d_calc_final["opt_constraint"]:
                    d["output"]["constraint"] = [
                        d_calc_final["opt_constraint"][0],
                        float(d_calc_final["opt_constraint"][6]),
                    ]
            if d["output"]["job_type"] in ["freq", "frequency"]:
                d["output"]["frequencies"] = d_calc_final["frequencies"]
                d["output"]["frequency_modes"] = d_calc_final["frequency_mode_vectors"]
                d["output"]["enthalpy"] = d_calc_final["total_enthalpy"]
                d["output"]["entropy"] = d_calc_final["total_entropy"]
                if d["input"]["job_type"] in ["opt", "optimization", "ts"]:
                    d["output"]["optimized_molecule"] = d_calc_final["initial_molecule"]
                    d["output"]["final_energy"] = d["calcs_reversed"][1]["final_energy"]
                # For frequency-first calcs
                else:
                    if len(d["calcs_reversed"]) > 1:
                        first_calc = d["calcs_reversed"][-1]["input"]["rem"]["job_type"]
                        second_calc = d["calcs_reversed"][-2]["input"]["rem"]["job_type"]
                        if first_calc in ["freq", "frequency"] and second_calc in ["opt", "optimization", "ts"]:
                            d["output"]["optimized_molecule"] = d_calc_final["initial_molecule"]
                            d["output"]["final_energy"] = d["calcs_reversed"][1]["final_energy"]

            if d["output"]["job_type"] == "pes_scan":
                d["output"]["scan_energies"] = d_calc_final.get("scan_energies")
                d["input"]["scan_variables"] = d_calc_final.get("scan_variables")
                d["output"]["scan_geometries"] = d_calc_final.get("optimized_geometries")
                d["output"]["scan_molecules"] = d_calc_final.get("molecules_from_optimized_geometries")

            if d["output"]["job_type"] == "force":
                d["output"]["gradients"] = d_calc_final["gradients"][0]
                if d_calc_final["pcm_gradients"] is not None:
                    d["output"]["pcm_gradients"] = d_calc_final["pcm_gradients"][0]
                if d_calc_final["CDS_gradients"] is not None:
                    d["output"]["CDS_gradients"] = d_calc_final["CDS_gradients"][0]

            opt_trajectory = []
            calcs = copy.deepcopy(d["calcs_reversed"])
            calcs.reverse()
            for calc in calcs:
                job_type = calc["input"]["rem"]["job_type"]
                if job_type in ["opt", "optimization", "ts"]:
                    for ii,geom in enumerate(calc["geometries"]):
                        site_properties = {"Mulliken":calc["Mulliken"][ii]}
                        if "RESP" in calc:
                            site_properties["RESP"] = calc["RESP"][ii]
                        mol = Molecule(
                            species=calc["species"],
                            coords=geom,
                            charge=calc["charge"],
                            spin_multiplicity=calc["multiplicity"],
                            site_properties=site_properties,
                        )
                        traj_entry = {"molecule": mol}
                        traj_entry["energy"] = calc["energy_trajectory"][ii]
                        opt_trajectory.append(traj_entry)
            if opt_trajectory != []:
                d["opt_trajectory"] = opt_trajectory

            if "final_energy" not in d["output"]:
                if d_calc_final["final_energy"] is not None:
                    d["output"]["final_energy"] = d_calc_final["final_energy"]
                else:
                    d["output"]["final_energy"] = d_calc_final["SCF"][-1][-1][0]

            if d_calc_final["completion"]:
                total_cputime = 0.0
                total_walltime = 0.0
                for calc in d["calcs_reversed"]:
                    if calc["walltime"] is not None:
                        total_walltime += calc["walltime"]
                    if calc["cputime"] is not None:
                        total_cputime += calc["cputime"]
                d["walltime"] = total_walltime
                d["cputime"] = total_cputime
            else:
                d["walltime"] = None
                d["cputime"] = None

            comp = d["output"]["initial_molecule"].composition
            d["formula_pretty"] = comp.reduced_formula
            d["formula_anonymous"] = comp.anonymized_formula
            d["formula_alphabetical"] = comp.alphabetical_formula
            d["chemsys"] = "-".join(sorted(set(d_calc_final["species"])))
            if d_calc_final["point_group"] is not None:
                d["pointgroup"] = d_calc_final["point_group"]
            else:
                try:
                    d["pointgroup"] = PointGroupAnalyzer(
                        d["output"]["initial_molecule"]
                    ).sch_symbol
                except ValueError:
                    d["pointgroup"] = "PGA_error"

            bb = BabelMolAdaptor(d["output"]["initial_molecule"])
            pbmol = bb.pybel_mol
            smiles = pbmol.write("smi").split()[0]
            d["smiles"] = smiles

            d["state"] = "successful" if d_calc_final["completion"] else "unsuccessful"
            if "special_run_type" in d:
                if d["special_run_type"] in ["frequency_flattener", "ts_frequency_flattener"]:
                    if d["state"] == "successful":
                        orig_num_neg_freq = None
                        for calc in d["calcs_reversed"][::-1]:
                            if "frequencies" in calc:
                                orig_num_neg_freq = sum(1 for freq in calc["frequencies"] if freq < 0)
                                break
                        orig_energy = None
                        for calc in d["calcs_reversed"][::-1]:
                            if "final_energy" in calc:
                                if calc["final_energy"] is not None:
                                    orig_energy = calc["final_energy"]
                                    break
                        final_num_neg_freq = sum(1 for freq in d_calc_final["frequencies"] if freq < 0)
                        final_energy = d["calcs_reversed"][1]["final_energy"]
                        if orig_num_neg_freq is None:
                            d["num_frequencies_flattened"] = 0
                        else:
                            d["num_frequencies_flattened"] = orig_num_neg_freq - final_num_neg_freq

                        if d["special_run_type"] == "frequency_flattener":
                            if final_num_neg_freq > 0: # If a negative frequency remains,
                                # and it's too large to ignore,
                                if final_num_neg_freq > 1 or abs(d["output"]["frequencies"][0]) >= 15.0:
                                    d["state"] = "unsuccessful" # then the flattening was unsuccessful
                        else:
                            if final_num_neg_freq > 1: # If a negative frequency remains,
                                # and it's too large to ignore,
                                if final_num_neg_freq > 2 or abs(d["output"]["frequencies"][1]) >= 15.0:
                                    d["state"] = "unsuccessful" # then the flattening was unsuccessful
                        if final_energy > orig_energy:
                            d["warnings"]["energy_increased"] = True

            d["last_updated"] = datetime.datetime.utcnow()
            return d

        except Exception:
            logger.error(traceback.format_exc())
            logger.error(
                "Error in " + os.path.abspath(dir_name) + ".\n" + traceback.format_exc()
            )
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
        d["input"]["smx"] = temp_input.smx
        d["input"]["vdw_mode"] = temp_input.vdw_mode
        d["input"]["van_der_waals"] = temp_input.van_der_waals
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
                "ERROR: The drone can only process a directory containing a single input/output pair when each "
                "include multiple calculations."
            )
        else:
            for key in input_files:
                to_return = []
                qchem_input_file = os.path.join(dir_name, input_files.get(key))
                qchem_output_file = os.path.join(dir_name, output_files.get(key))
                multi_out = QCOutput.multiple_outputs_from_file(
                    QCOutput, qchem_output_file, keep_sub_files=False
                )
                multi_in = QCInput.from_multi_jobs_file(qchem_input_file)
                for ii, out in enumerate(multi_out):
                    d = out.data
                    d["input"] = {}
                    d["input"]["molecule"] = multi_in[ii].molecule
                    d["input"]["rem"] = multi_in[ii].rem
                    d["input"]["opt"] = multi_in[ii].opt
                    d["input"]["pcm"] = multi_in[ii].pcm
                    d["input"]["solvent"] = multi_in[ii].solvent
                    d["input"]["smx"] = multi_in[ii].smx
                    d["input"]["vdw_mode"] = multi_in[ii].vdw_mode
                    d["input"]["van_der_waals"] = multi_in[ii].van_der_waals
                    d["task"] = {"type": key, "name": "calc" + str(ii)}
                    to_return.append(d)
            return to_return

    @staticmethod
    def post_process(dir_name, d):
        """
        Post-processing for various files other than the QChem input and output files.
        """
        logger.info(f"Post-processing dir:{dir_name}")
        fullpath = os.path.abspath(dir_name)
        filenames = glob.glob(os.path.join(fullpath, "custodian.json*"))
        if len(filenames) >= 1:
            with zopen(filenames[0], "rt") as f:
                d["custodian"] = json.load(f)
        filenames = glob.glob(os.path.join(fullpath, "solvent_data*"))
        if len(filenames) >= 1:
            with zopen(filenames[0], "rt") as f:
                d["custom_smd"] = f.readlines()[0]
        filenames = glob.glob(os.path.join(fullpath, "processed_critic2.json*"))
        if len(filenames) >= 1:
            with zopen(filenames[0], "rt") as f:
                d["critic2"] = {}
                d["critic2"]["processed"] = json.load(f)
            filenames = glob.glob(os.path.join(fullpath, "cpreport.json*"))
            if len(filenames) >= 1:
                with zopen(filenames[0], "rt") as f:
                    d["critic2"]["CP"] = json.load(f)
            filenames = glob.glob(os.path.join(fullpath, "yt.json*"))
            if len(filenames) >= 1:
                with zopen(filenames[0], "rt") as f:
                    d["critic2"]["YT"] = json.load(f)
            filenames = glob.glob(os.path.join(fullpath, "bonding.json*"))
            if len(filenames) >= 1:
                with zopen(filenames[0], "rt") as f:
                    d["critic2"]["bonding"] = json.load(f)

    def validate_doc(self, d):
        """
        Sanity check, aka make sure all the important keys are set. Note that a failure
        to pass validation is unfortunately unlikely to be noticed by a user.
        """
        for k, v in self.schema.items():
            diff = v.difference(set(d.get(k, d).keys()))
            if diff:
                logger.warning(f"The keys {diff} in {k} not set")

    @staticmethod
    def get_valid_paths(self, path):
        return [path]
