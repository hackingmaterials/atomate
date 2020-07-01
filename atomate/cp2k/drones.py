# coding: utf-8

"""
This drone assimilates directories containing the results of cp2k calculations
"""

import os
import re
import datetime
from fnmatch import fnmatch
from collections import OrderedDict
import json
import glob
import traceback
import warnings

from monty.io import zopen
from monty.json import jsanitize

import numpy as np

from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.cp2k.outputs import Cp2kOutput, Cube
from pymatgen.apps.borg.hive import AbstractDrone
from pymatgen.io.vasp.outputs import VolumetricData

from atomate.utils.utils import get_logger, get_uri
from atomate import __version__ as atomate_version

__author__ = "Nicholas Winner"

logger = get_logger(__name__)


# TODO: Don't push until going through all this
class Cp2kDrone(AbstractDrone):
    """
    Adapted from VaspDrone
    """

    __version__ = (
        atomate_version  # note: the version is inserted into the task doc
    )

    # Schema def of important keys and sub-keys; used in validation
    schema = {
        "root": {
            "schema",
            "dir_name",
            "chemsys",
            "composition_reduced",
            "formula_pretty",
            "formula_reduced_abc",
            "elements",
            "nelements",
            "formula_anonymous",
            "calcs_reversed",
            "completed_at",
            "nsites",
            "composition_unit_cell",
            "input",
            "output",
            "state",
            "analysis",
            "run_stats",
        },
        "input": {"atomic_kind_info", "structure"},
        "output": {
            "structure",
            "spacegroup",
            "density",
            "energy",
            "energy_per_atom",
            "bandgap",
            "vbm",
            "cbm",
            "is_metal",
            "forces",
            "stress",
            "ionic_steps"
        },
        "calcs_reversed": {
            "dir_name",
            "run_type",
            "elements",
            "nelements",
            "formula_pretty",
            "formula_reduced_abc",
            "composition_reduced",
            "cp2k_version",
            "formula_anonymous",
            "nsites",
            "composition_unit_cell",
            "completed_at",
            "task",
            "input",
            "output",
            "has_cp2k_completed",
        },
        "analysis": {
            "delta_volume_as_percent",
            "delta_volume",
            "max_force",
            "errors",
            "warnings",
        },
    }

    def __init__(
        self,
        runs=None,
        parse_dos="auto",
        parse_hartree=False,
        additional_fields=None,
        use_full_uri=True,
    ):
        """
        Initialize a cp2k drone to parse cp2k outputs
        Args:
            runs (list): Naming scheme for multiple calcuations in one folder e.g. ["relax1","relax2"].
             Can be subfolder or extension
            parse_dos (str or bool): Whether to parse the DOS. Can be "auto", True or False.
                "auto" will only parse DOS if NSW = 0, so there are no ionic steps

            additional_fields (dict): dictionary of additional fields to add to output document
            use_full_uri (bool): converts the directory path to the full URI path
        """
        self.runs = runs or ["precondition"] + [
            "relax" + str(i + 1) for i in range(9)
        ]
        self.parse_dos = parse_dos
        self.parse_hartree = parse_hartree
        self.additional_fields = additional_fields or {}
        self.use_full_uri = use_full_uri

    def assimilate(self, path):
        """
        Adapted from matgendb.creator
        Parses cp2k runs and insert the result into the db.
        Get the entire task doc from the cp2k.out file in the path.
        Also adds some post-processed info.

        Args:
            path (str): Path to the directory containing cp2k.out file

        Returns:
            (dict): a task dictionary
        """
        logger.info("Getting task doc for base dir :{}".format(path))
        cp2k_out_files = self.filter_files(path, file_pattern="cp2k.out")
        if len(cp2k_out_files) > 0:
            d = self.generate_doc(path, cp2k_out_files)
            self.post_process(path, d)
        else:
            raise ValueError("No CP2K output files found!")
        self.validate_doc(d)
        return d

    def filter_files(self, path, file_pattern="cp2k.out"):
        """
        Find the files that match the pattern in the given path and
        return them in an ordered dictionary. The searched for files are
        filtered by the run types defined in self.runs. e.g. ["relax1", "relax2", ...].
        Only 2 schemes of the file filtering is enabled: searching for run types
        in the list of files and in the filenames. Modify this method if more
        sophisticated filtering scheme is needed.

        Args:
            path (string): path to the folder
            file_pattern (string): files to be searched for

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

    def generate_doc(self, dir_name, cp2k_files):
        """
        Adapted from matgendb.creator.generate_doc
        """
        try:
            # basic properties, incl. calcs_reversed and run_stats
            fullpath = os.path.abspath(dir_name)
            d = jsanitize(self.additional_fields, strict=True)
            d["schema"] = {"code": "atomate", "version": Cp2kDrone.__version__}
            d["dir_name"] = fullpath
            d["calcs_reversed"] = [
                self.process_cp2k(dir_name, taskname, filename)
                for taskname, filename in cp2k_files.items()
            ]

            # TODO More run stats
            run_stats = {}
            for i, d_calc in enumerate(d["calcs_reversed"]):
                run_stats[d_calc["task"]["name"]] = d_calc["total_time"]
            d["run_stats"] = run_stats

            # TODO
            # set root formula/composition keys based on initial and final calcs
            d_calc_init = d["calcs_reversed"][-1]
            d_calc_final = d["calcs_reversed"][0]
            d["chemsys"] = "-".join(sorted(d_calc_final["elements"]))
            comp = Composition(d_calc_final["composition_unit_cell"])
            d["formula_anonymous"] = comp.anonymized_formula
            d[
                "formula_reduced_abc"
            ] = comp.reduced_composition.alphabetical_formula
            for root_key in [
                "completed_at",
                "nsites",
                "composition_unit_cell",
                "composition_reduced",
                "formula_pretty",
                "elements",
                "nelements",
            ]:
                d[root_key] = d_calc_final[root_key]

            # TODO: Should atomic kind info and DFT really be saved like this?
            # Right now the input is where you save the dft parameters, but for
            # things like hybrid calcluations, the dft parameters change as you
            # switch from GGA to GGA hybrid. Is it maybe better to save dft
            # parameters as a sequence or something? showing how they are the whole time?
            d["input"] = {
                "structure": d_calc_init["input"]["structure"],
                "atomic_kind_info": d_calc_init["input"]["atomic_kind_info"],
                "dft": d_calc_init["input"]["dft"],
                "scf": d_calc_init["input"]["scf"],
            }

            # store the output key based on final calc
            d["output"] = {
                "structure": d_calc_final["output"]["structure"],
                "density": d_calc_final.pop("density"),
                "energy": d_calc_final["output"]["energy"],
                "energy_per_atom": d_calc_final["output"]["energy_per_atom"],
                "forces": d_calc_final["output"]["ionic_steps"][-1].get(
                    "forces"
                ),
                "stress": d_calc_final["output"]["ionic_steps"][-1].get(
                    "stress"
                ),
                "ionic_steps": d_calc_final["output"]["ionic_steps"],
                "cbm": d_calc_final["output"]["cbm"],
                "vbm": d_calc_final["output"]["vbm"],
                "bandgap": d_calc_final["output"]["bandgap"],
                "efermi": d_calc_final["output"]["efermi"],
                "is_metal": d_calc_final["output"]["is_metal"],
                "v_hartree": d_calc_final.pop('v_hartree', None),
                "v_hartree_grid": d_calc_final.pop('v_hartree_grid', None),
            }

            # Store symmetry information
            try:
                sg = SpacegroupAnalyzer(
                    Structure.from_dict(d_calc_final["output"]["structure"]),
                    0.1,
                )
                if not sg.get_symmetry_dataset():
                    sg = SpacegroupAnalyzer(
                        Structure.from_dict(
                            d_calc_final["output"]["structure"]
                        ),
                        1e-3,
                        1,
                    )
                d["output"]["spacegroup"] = {
                    "source": "spglib",
                    "symbol": sg.get_space_group_symbol(),
                    "number": sg.get_space_group_number(),
                    "point_group": sg.get_point_group_symbol(),
                    "crystal_system": sg.get_crystal_system(),
                    "hall": sg.get_hall(),
                }
            except TypeError:
                d["output"]["spacegroup"] = {
                    "source": None,
                    "symbol": None,
                    "number": None,
                    "point_group": None,
                    "crystal_system": None,
                    "hall": None,
                }
                warnings.warn(
                    "Space Group could not be determined by this drone.",
                    Warning,
                )

            d["state"] = (
                "successful"
                if all([i["has_cp2k_completed"] for i in d["calcs_reversed"]])
                else "unsuccessful"
            )

            self.set_analysis(d)

            d["last_updated"] = datetime.datetime.utcnow()
            return d

        except Exception:
            logger.error(traceback.format_exc())
            logger.error(
                "Error in "
                + os.path.abspath(dir_name)
                + ".\n"
                + traceback.format_exc()
            )
            raise

    def process_cp2k(self, dir_name, taskname, filename):
        """
        Adapted from matgendb.creator

        Process a cp2k output file.
        """
        cp2k_file = os.path.join(dir_name, filename)

        out = Cp2kOutput(cp2k_file, auto_load=True)
        d = out.as_dict()

        comp = Composition(d["composition"])
        d["formula_pretty"] = comp.reduced_formula
        d["composition_reduced"] = comp.reduced_composition
        d["composition_unit_cell"] = comp.as_dict()
        d["formula_anonymous"] = comp.anonymized_formula
        d["formula_reduced_abc"] = comp.reduced_composition.alphabetical_formula
        d["elements"] = list(comp.as_dict().keys())
        d["nelements"] = len(self.as_dict().keys())
        d["nsites"] = len(d["input"]["structure"]["sites"])
        d["dir_name"] = os.path.abspath(dir_name)
        d["completed_at"] = str(
            datetime.datetime.fromtimestamp(os.path.getmtime(cp2k_file))
        )
        d["density"] = out.final_structure.density

        d["has_cp2k_completed"] = d.pop("ran_successfully")

        # store run name and location ,e.g. relax1, relax2, etc.
        d["task"] = {"type": taskname, "name": taskname}

        # include output file names
        d["output_file_paths"] = self.process_raw_data(
            dir_name, taskname=taskname
        )

        if self.parse_dos:
            d['dos'] = out.parse_pdos()

        if self.parse_hartree:
            cube = Cube(out.filenames['v_hartree'][-1])
            vd = VolumetricData(structure=cube.structure, data={'total': cube.data})
            d['v_hartree'] = [
                    vd.get_average_along_axis(i) for i in range(3)
            ]
            d['v_hartree_grid'] = [
                vd.get_axis_grid(i) for i in range(3)
            ]
        return d

    def process_raw_data(self, dir_name, taskname="standard"):
        pass

    @staticmethod
    def set_analysis(d, max_force_threshold=0.5, volume_change_threshold=0.2):
        """
        Adapted from matgendb.creator

        set the 'analysis' key
        """
        initial_vol = d["input"]["structure"]["lattice"]["volume"]
        final_vol = d["output"]["structure"]["lattice"]["volume"]
        delta_vol = final_vol - initial_vol
        percent_delta_vol = 100 * delta_vol / initial_vol
        warning_msgs = []
        error_msgs = []

        # delta volume checks
        if abs(percent_delta_vol) > volume_change_threshold:
            warning_msgs.append(
                "Volume change > {}%".format(volume_change_threshold * 100)
            )

        # max force and valid structure checks
        max_force = None
        calc = d["calcs_reversed"][0]
        if (
            d["state"] == "successful"
            and "forces" in calc["output"]["ionic_steps"][-1].keys()
        ):

            # calculate max forces
            forces = np.array(calc["output"]["ionic_steps"][-1]["forces"])
            # account for selective dynamics
            final_structure = Structure.from_dict(calc["output"]["structure"])
            sdyn = final_structure.site_properties.get("selective_dynamics")
            if sdyn:
                forces[np.logical_not(sdyn)] = 0
            max_force = max(np.linalg.norm(forces, axis=1))

            s = Structure.from_dict(d["output"]["structure"])
            if not s.is_valid():
                error_msgs.append("Bad structure (atoms are too close!)")
                d["state"] = "error"

        d["analysis"] = {
            "delta_volume": delta_vol,
            "delta_volume_as_percent": percent_delta_vol,
            "max_force": max_force,
            "warnings": warning_msgs,
            "errors": error_msgs,
        }

    def post_process(self, dir_name, d):
        """
        Post-processing for various files other than the cp2k.out file.
        Looks for files: transformations.json and custodian.json. Modify this if other
        output files need to be processed.

        Args:
            dir_name:
                The dir_name.
            d:
                Current doc generated.
        """
        logger.info("Post-processing dir:{}".format(dir_name))
        fullpath = os.path.abspath(dir_name)
        # CP2K input generated by pymatgen's alchemy has a transformations.json file that tracks
        # the origin of a particular structure. If such a file is found, it is inserted into the
        # task doc as d["transformations"]
        transformations = {}
        filenames = glob.glob(os.path.join(fullpath, "transformations.json*"))
        if len(filenames) >= 1:
            with zopen(filenames[0], "rt") as f:
                transformations = json.load(f)
                try:
                    m = re.match(
                        "(\d+)-ICSD", transformations["history"][0]["source"]
                    )
                    if m:
                        d["icsd_id"] = int(m.group(1))
                except Exception as ex:
                    logger.warning(
                        "Cannot parse ICSD from transformations file."
                    )
                    pass
        else:
            logger.warning("Transformations file does not exist.")

        other_parameters = transformations.get("other_parameters")
        new_tags = None
        if other_parameters:
            # We don't want to leave tags or authors in the
            # transformations file because they'd be copied into
            # every structure generated after this one.
            new_tags = other_parameters.pop("tags", None)
            new_author = other_parameters.pop("author", None)
            if new_author:
                d["author"] = new_author
            if not other_parameters:  # if dict is now empty remove it
                transformations.pop("other_parameters")
        d["transformations"] = transformations

        # Calculations done using custodian has a custodian.json,
        # which tracks the jobs performed and any errors detected and fixed.
        # This is useful for tracking what has actually be done to get a
        # result. If such a file is found, it is inserted into the task doc
        # as d["custodian"]
        filenames = glob.glob(os.path.join(fullpath, "custodian.json*"))
        if len(filenames) >= 1:
            with zopen(filenames[0], "rt") as f:
                d["custodian"] = json.load(f)
        # Convert to full uri path.
        if self.use_full_uri:
            d["dir_name"] = get_uri(dir_name)
        if new_tags:
            d["tags"] = new_tags

        # Calculations using custodian generate a *.orig file for the inputs
        # This is useful to know how the calculation originally started
        # if such files are found they are inserted into orig_inputs
        filenames = glob.glob(os.path.join(fullpath, "*.orig*"))

        if len(filenames) >= 1:
            d["orig_inputs"] = {}
            for f in filenames:
                if "INCAR.orig" in f:
                    d["orig_inputs"]["incar"] = Incar.from_file(f).as_dict()
                if "POTCAR.orig" in f:
                    d["orig_inputs"]["potcar"] = Potcar.from_file(f).as_dict()
                if "KPOINTS.orig" in f:
                    d["orig_inputs"]["kpoints"] = Kpoints.from_file(f).as_dict()
                if "POSCAR.orig" in f:
                    d["orig_inputs"]["poscar"] = Poscar.from_file(f).as_dict()

        logger.info("Post-processed " + fullpath)

    def validate_doc(self, d):
        """
        Sanity check.
        Make sure all the important keys are set
        """
        # TODO: @matk86 - I like the validation but I think no one will notice a failed
        # validation tests which removes the usefulness of this. Any ideas to make people
        # notice if the validation fails? -computron
        for k, v in self.schema.items():
            if k == "calcs_reversed":
                diff = v.difference(set(d.get(k, d)[0].keys()))
            else:
                diff = v.difference(set(d.get(k, d).keys()))
            if diff:
                logger.warn("The keys {0} in {1} not set".format(diff, k))

    def get_valid_paths(self, path):
        """
        There are some restrictions on the valid directory structures:

        1. There can be only one cp2k run in each directory. Nested directories
           are fine.
        2. Directories designated "relax1"..."relax9" are considered to be
           parts of a multiple-optimization run.
        3. Directories containing cp2k output with ".relax1"...".relax9" are
           also considered as parts of a multiple-optimization run.
        """
        (parent, subdirs, files) = path
        if set(self.runs).intersection(subdirs):
            return [parent]
        if (
            not any([parent.endswith(os.sep + r) for r in self.runs])
            and len(glob.glob(os.path.join(parent, "cp2k.out*"))) > 0
        ):
            return [parent]
        return []

    def as_dict(self):
        init_args = {
            "parse_dos": self.parse_dos,
            "additional_fields": self.additional_fields,
            "use_full_uri": self.use_full_uri,
            "runs": self.runs,
        }
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "version": self.__class__.__version__,
            "init_args": init_args,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(**d["init_args"])
