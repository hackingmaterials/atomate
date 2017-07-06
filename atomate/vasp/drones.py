# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This Drone tries to produce a more sensible task dictionary than the default VaspToDbTaskDrone.
Some of the changes are documented in this thread:
https://groups.google.com/forum/#!topic/pymatgen/pQ-emBpeV5U
"""

import os
import re
import datetime
from fnmatch import fnmatch
from collections import OrderedDict
import json
import glob

from monty.io import zopen

import numpy as np

from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Vasprun, Outcar
from pymatgen.apps.borg.hive import AbstractDrone

from matgendb.creator import get_uri

from atomate.utils.utils import get_logger

__author__ = 'Kiran Mathew, Shyue Ping Ong, Shyam Dwaraknath, Anubhav Jain'
__email__ = 'kmathew@lbl.gov'
__date__ = 'Mar 27, 2016'
__version__ = "0.1.0"

logger = get_logger(__name__)


class VaspDrone(AbstractDrone):
    """
    pymatgen-db VaspToDbTaskDrone with updated schema and documents processing methods.
    Please refer to matgendb.creator.VaspToDbTaskDrone documentation.
    """

    __version__ = 0.2  # note: the version is inserted into the task doc

    # Schema def of important keys and sub-keys; used in validation
    schema = {
        "root": {
            "schema", "dir_name", "chemsys", "composition_reduced",
            "formula_pretty", "formula_reduced_abc", "elements",
            "nelements", "formula_anonymous", "calcs_reversed", "completed_at",
            "nsites", "composition_unit_cell", "input", "output", "state",
            "analysis", "run_stats"
        },
        "input": {'is_lasph', 'is_hubbard', 'xc_override', 'potcar_spec',
                  'hubbards', 'structure', 'pseudo_potential'},
        "output": {'structure', 'spacegroup', 'density', 'energy',
                   'energy_per_atom', 'is_gap_direct', 'bandgap', 'vbm',
                   'cbm', 'is_metal'},
        "calcs_reversed": {
            'dir_name', 'run_type', 'elements', 'nelements',
            'formula_pretty', 'formula_reduced_abc', 'composition_reduced',
            'vasp_version', 'formula_anonymous', 'nsites',
            'composition_unit_cell', 'completed_at', 'task', 'input', 'output',
            'has_vasp_completed'
        },
        "analysis": {'delta_volume_percent', 'delta_volume', 'max_force',
                     'errors',
                     'warnings'}
    }

    def __init__(self, runs=None, parse_dos=False, compress_dos=False, bandstructure_mode=False,
                 compress_bs=False, additional_fields=None, use_full_uri=True):
        self.parse_dos = parse_dos
        self.compress_dos = compress_dos
        self.additional_fields = additional_fields or {}
        self.use_full_uri = use_full_uri
        self.runs = runs or ["relax" + str(i+1) for i in range(9)]  # can't auto-detect: path unknown
        self.bandstructure_mode = bandstructure_mode
        self.compress_bs = compress_bs

    def assimilate(self, path):
        """
        Adapted from matgendb.creator
        Parses vasp runs(vasprun.xml file) and insert the result into the db.
        Get the entire task doc from the vasprum.xml and the OUTCAR files in the path.
        Also adds some post-processed info.

        Args:
            path (str): Path to the directory containing vasprun.xml and OUTCAR files

        Returns:
            (dict): a task dictionary
        """
        logger.info("Getting task doc for base dir :{}".format(path))
        vasprun_files = self.filter_files(path, file_pattern="vasprun.xml")
        outcar_files = self.filter_files(path, file_pattern="OUTCAR")
        if len(vasprun_files) > 0 and len(outcar_files) > 0:
            d = self.generate_doc(path, vasprun_files, outcar_files)
            self.post_process(path, d)
        else:
            raise ValueError("No VASP files found!")
        self.validate_doc(d)
        return d

    def filter_files(self, path, file_pattern="vasprun.xml"):
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
                    processed_files['standard'] = f
        return processed_files

    def generate_doc(self, dir_name, vasprun_files, outcar_files):
        """
        Adapted from matgendb.creator.generate_doc
        """
        try:
            # basic properties, incl. calcs_reversed and run_stats
            fullpath = os.path.abspath(dir_name)
            d = {k: v for k, v in self.additional_fields.items()}
            d["schema"] = {"code": "atomate", "version": VaspDrone.__version__}
            d["dir_name"] = fullpath
            d["calcs_reversed"] = [self.process_vasprun(dir_name, taskname, filename)
                                   for taskname, filename in vasprun_files.items()]
            outcar_data = [Outcar(os.path.join(dir_name, filename)).as_dict()
                           for taskname, filename in outcar_files.items()]
            run_stats = {}
            for i, d_calc in enumerate(d["calcs_reversed"]):
                run_stats[d_calc["task"]["name"]] = outcar_data[i].pop("run_stats")
                if d_calc.get("output"):
                    d_calc["output"].update({"outcar": outcar_data[i]})
                else:
                    d_calc["output"] = {"outcar": outcar_data[i]}
            try:
                overall_run_stats = {}
                for key in ["Total CPU time used (sec)", "User time (sec)", "System time (sec)",
                            "Elapsed time (sec)"]:
                    overall_run_stats[key] = sum([v[key] for v in run_stats.values()])
                run_stats["overall"] = overall_run_stats
            except:
                logger.error("Bad run stats for {}.".format(fullpath))
            d["run_stats"] = run_stats

            # reverse the calculations data order so newest calc is first
            d["calcs_reversed"].reverse()

            # set root formula/composition keys based on initial and final calcs
            d_calc_init = d["calcs_reversed"][-1]
            d_calc_final = d["calcs_reversed"][0]
            d["chemsys"] = "-".join(sorted(d_calc_final["elements"]))
            comp = Composition(d_calc_final["composition_unit_cell"])
            d["formula_anonymous"] = comp.anonymized_formula
            d["formula_reduced_abc"] = comp.reduced_composition.alphabetical_formula
            for root_key in ["completed_at", "nsites", "composition_unit_cell",
                             "composition_reduced", "formula_pretty", "elements", "nelements"]:
                d[root_key] = d_calc_final[root_key]

            # store the input key based on initial calc
            # store any overrides to the exchange correlation functional
            xc = d_calc_init["input"]["incar"].get("GGA")
            if xc:
                xc = xc.upper()
            p = d_calc_init["input"]["potcar_type"][0].split("_")
            pot_type = p[0]
            functional = "lda" if len(pot_type) == 1 else "_".join(p[1:])
            d["input"] = {"structure": d_calc_init["input"]["structure"],
                          "is_hubbard": d_calc_init.pop("is_hubbard"),
                          "hubbards": d_calc_init.pop("hubbards"),
                          "is_lasph": d_calc_init["input"]["incar"].get("LASPH", False),
                          "potcar_spec": d_calc_init["input"].get("potcar_spec"),
                          "xc_override": xc,
                          "pseudo_potential": {"functional": functional.lower(),
                                               "pot_type": pot_type.lower(),
                                               "labels": d_calc_init["input"]["potcar"]},
                          "parameters": d_calc_init["input"]["parameters"],
                          "incar": d_calc_init["input"]["incar"]
                          }

            # store the output key based on final calc
            d["output"] = {
                "structure": d_calc_final["output"]["structure"],
                "density": d_calc_final.pop("density"),
                "energy": d_calc_final["output"]["energy"],
                "energy_per_atom": d_calc_final["output"]["energy_per_atom"]}
            
            # patch calculated magnetic moments into final structure
            if len(d_calc_final["output"]["outcar"]["magnetization"]) != 0:
                magmoms = [m["tot"] for m in d_calc_final["output"]["outcar"]["magnetization"]]
                s = Structure.from_dict(d["output"]["structure"])
                s.add_site_property('magmom', magmoms)
                d["output"]["structure"] = s.as_dict()

            calc = d["calcs_reversed"][0]
            gap = calc["output"]["bandgap"]
            cbm = calc["output"]["cbm"]
            vbm = calc["output"]["vbm"]
            is_direct = calc["output"]["is_gap_direct"]
            is_metal = calc["output"]["is_metal"]
            d["output"].update({"bandgap": gap, "cbm": cbm, "vbm": vbm,
                                "is_gap_direct": is_direct, "is_metal": is_metal})

            sg = SpacegroupAnalyzer(Structure.from_dict(d_calc_final["output"]["structure"]), 0.1)
            if not sg.get_symmetry_dataset():
                sg = SpacegroupAnalyzer(Structure.from_dict(d_calc_final["output"]["structure"]),
                                        1e-3, 1)
            d["output"]["spacegroup"] = {
                "source": "spglib",
                "symbol": sg.get_space_group_symbol(),
                "number": sg.get_space_group_number(),
                "point_group": sg.get_point_group_symbol(),
                "crystal_system": sg.get_crystal_system(),
                "hall": sg.get_hall()}
            if d["input"]["parameters"].get("LEPSILON"):
                for k in ['epsilon_static', 'epsilon_static_wolfe', 'epsilon_ionic']:
                    d["output"][k] = d_calc_final["output"][k]

            d["state"] = "successful" if d_calc["has_vasp_completed"] else "unsuccessful"

            self.set_analysis(d)

            d["last_updated"] = datetime.datetime.today()
            return d

        except Exception:
            import traceback
            logger.error(traceback.format_exc())
            logger.error("Error in " + os.path.abspath(dir_name) + ".\n" + traceback.format_exc())
            raise

    def process_vasprun(self, dir_name, taskname, filename):
        """
        Adapted from matgendb.creator

        Process a vasprun.xml file.
        """
        vasprun_file = os.path.join(dir_name, filename)
        if self.bandstructure_mode:
            vrun = Vasprun(vasprun_file, parse_eigen=True, parse_projected_eigen=True)
        else:
            vrun = Vasprun(vasprun_file)

        d = vrun.as_dict()

        # rename formula keys
        for k, v in {"formula_pretty": "pretty_formula",
                     "composition_reduced": "reduced_cell_formula",
                     "composition_unit_cell": "unit_cell_formula"}.items():
            d[k] = d.pop(v)

        for k in ["eigenvalues", "projected_eigenvalues"]:  # large storage space breaks some docs
            if k in d["output"]:
                del d["output"][k]

        comp = Composition(d["composition_unit_cell"])
        d["formula_anonymous"] = comp.anonymized_formula
        d["formula_reduced_abc"] = comp.reduced_composition.alphabetical_formula
        d["dir_name"] = os.path.abspath(dir_name)
        d["completed_at"] = str(datetime.datetime.fromtimestamp(os.path.getmtime(vasprun_file)))
        d["density"] = vrun.final_structure.density

        # replace 'crystal' with 'structure'
        d["input"]["structure"] = d["input"].pop("crystal")
        d["output"]["structure"] = d["output"].pop("crystal")
        for k, v in {"energy": "final_energy", "energy_per_atom": "final_energy_per_atom"}.items():
            d["output"][k] = d["output"].pop(v)

        if self.parse_dos and self.parse_dos != 'final':
            try:
                d["dos"] = vrun.complete_dos.as_dict()
            except:
                raise ValueError("No valid dos data exist in {}.".format(dir_name))

        if self.bandstructure_mode:
            bs = vrun.get_band_structure(line_mode=(self.bandstructure_mode.lower() == "line"))
        else:
            bs = vrun.get_band_structure()

        d["bandstructure"] = bs.as_dict()

        d["output"]["vbm"] = bs.get_vbm()["energy"]
        d["output"]["cbm"] = bs.get_cbm()["energy"]
        bs_gap = bs.get_band_gap()
        d["output"]["bandgap"] = bs_gap["energy"]
        d["output"]["is_gap_direct"] = bs_gap["direct"]
        d["output"]["is_metal"] = bs.is_metal()
        d["task"] = {"type": taskname, "name": taskname}

        if hasattr(vrun, "force_constants"):
            # phonon-dfpt
            d["output"]["force_constants"] = vrun.force_constants.tolist()
            d["output"]["normalmode_eigenvals"] = vrun.normalmode_eigenvals.tolist()
            d["output"]["normalmode_eigenvecs"] = vrun.normalmode_eigenvecs.tolist()
        return d

    @staticmethod
    def set_analysis(d, max_force_threshold=0.5, volume_change_threshold=0.2):
        """
        Adapted from matgendb.creator

        set the 'analysis' key
        """
        initial_vol = d["input"]["structure"]["lattice"]["volume"]
        final_vol = d["output"]["structure"]["lattice"]["volume"]
        delta_vol = final_vol - initial_vol
        percent_delta_vol = delta_vol / initial_vol
        warning_msgs = []
        error_msgs = []

        # delta volume checks
        if abs(percent_delta_vol) > volume_change_threshold:
            warning_msgs.append("Volume change > {}%".format(volume_change_threshold * 100))

        # max force and valid structure checks
        max_force = None
        calc = d["calcs_reversed"][0]
        if d["state"] == "successful" and calc["input"]["parameters"].get("NSW", 0) > 0:
            # handle the max force and max force error
            max_force = max([np.linalg.norm(a) for a in calc["output"]["ionic_steps"][-1]["forces"]])
            if max_force > max_force_threshold:
                error_msgs.append("Final max force exceeds {} eV".format(max_force_threshold))
                d["state"] = "error"

            s = Structure.from_dict(d["output"]["structure"])
            if not s.is_valid():
                error_msgs.append("Bad structure (atoms are too close!)")
                d["state"] = "error"

        d["analysis"] = {"delta_volume": delta_vol,
                         "delta_volume_percent": percent_delta_vol,
                         "max_force": max_force,
                         "warnings": warning_msgs,
                         "errors": error_msgs}

    def post_process(self, dir_name, d):
        """
        Post-processing for various files other than the vasprun.xml and OUTCAR.
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
        # VASP input generated by pymatgen's alchemy has a transformations.json file that tracks
        # the origin of a particular structure. If such a file is found, it is inserted into the
        # task doc as d["transformations"]
        transformations = {}
        filenames = glob.glob(os.path.join(fullpath, "transformations.json*"))
        if len(filenames) >= 1:
            with zopen(filenames[0], "rt") as f:
                transformations = json.load(f)
                try:
                    m = re.match("(\d+)-ICSD", transformations["history"][0]["source"])
                    if m:
                        d["icsd_id"] = int(m.group(1))
                except Exception as ex:
                    logger.warning("Cannot parse ICSD from transformations file.")
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

        1. There can be only one vasp run in each directory. Nested directories
           are fine.
        2. Directories designated "relax1"..."relax9" are considered to be
           parts of a multiple-optimization run.
        3. Directories containing vasp output with ".relax1"...".relax9" are
           also considered as parts of a multiple-optimization run.
        """
        (parent, subdirs, files) = path
        if set(self.runs).intersection(subdirs):
            return [parent]
        if not any([parent.endswith(os.sep + r) for r in self.runs]) and \
                len(glob.glob(os.path.join(parent, "vasprun.xml*"))) > 0:
            return [parent]
        return []

    def as_dict(self):
        init_args = {
            "parse_dos": self.parse_dos,
            "compress_dos": self.compress_dos,
            "bandstructure_mode": self.bandstructure_mode,
            "compress_bs": self.compress_bs,
            "additional_fields": self.additional_fields,
            "use_full_uri": self.use_full_uri,
            "runs": self.runs}
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "version": self.__class__.__version__,
                "init_args": init_args
                }

    @classmethod
    def from_dict(cls, d):
        return cls(**d["init_args"])
