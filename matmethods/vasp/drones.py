# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This Drone tries to produce a more sensible task dictionary than
the default VaspToDbTaskDrone. Some of the changes are documented
in this thread:
https://groups.google.com/forum/#!topic/pymatgen/pQ-emBpeV5U
"""

import os
import re
import datetime
import zlib
from fnmatch import fnmatch
from collections import OrderedDict
import json
import glob

from monty.io import zopen
from monty.json import MontyEncoder

import numpy as np

from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Vasprun, Outcar
from pymatgen.apps.borg.hive import AbstractDrone
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.matproj.rest import MPRester

from pymongo import MongoClient
import gridfs

from matgendb.creator import get_uri

from matmethods.utils.utils import get_logger

__author__ = 'Kiran Mathew, Shyue Ping Ong, Shyam Dwaraknath'
__credits__ = 'Anubhav Jain'
__email__ = 'kmathew@lbl.gov'
__date__ = 'Mar 27, 2016'
__version__ = "0.1.0"

logger = get_logger(__name__)

# TODO: needs comprehensive unit tests


class VaspToDbTaskDrone(AbstractDrone):
    """
    pymatgen-db VaspToDbTaskDrone with updated schema and documents processing methods.
    Please refer to matgendb.creator.VaspToDbTaskDrone documentation.
    """

    __version__ = 0.1

    def __init__(self, host="127.0.0.1", port=27017, database="vasp", collection="tasks",
                 user=None, password=None, simulate_mode=False, runs=["relax1", "relax2"],
                 parse_dos=False, compress_dos=False, bandstructure_mode=False,
                 compress_bs=False, additional_fields={}, update_duplicates=True,
                 mapi_key=None, use_full_uri=True):
        self.host = host
        self.database = database
        self.user = user
        self.password = password
        self.collection = collection
        self.port = port
        self.simulate = simulate_mode
        self.parse_dos = parse_dos
        self.compress_dos = compress_dos
        self.additional_fields = additional_fields
        self.update_duplicates = update_duplicates
        self.mapi_key = mapi_key
        self.use_full_uri = use_full_uri
        self.runs = runs
        self.bandstructure_mode = bandstructure_mode
        self.compress_bs = compress_bs
        if not simulate_mode:
            conn = MongoClient(self.host, self.port, j=True)
            db = conn[self.database]
            if self.user:
                db.authenticate(self.user, self.password)
            if db.counter.find({"_id": "taskid"}).count() == 0:
                db.counter.insert({"_id": "taskid", "c": 1})
        # set of important keys and sub-keys
        self.root_keys = {"schema", "dir_name", "chemsys", "composition_reduced",
                          "formula_pretty", "formula_reduced_abc", "elements", "nelements",
                          "formula_anonymous", "calcs_reversed", "completed_at", "nsites",
                          "composition_unit_cell", "input", "output", "state", "analysis",
                          "run_stats"}
        self.input_keys = {'is_lasph', 'is_hubbard', 'xc_override', 'potcar_spec', 'hubbards',
                           'structure', 'pseudo_potential'}
        self.output_keys = {'structure', 'spacegroup', 'density', 'energy', 'energy_per_atom',
                            'is_gap_direct',  'bandgap', 'vbm', 'cbm', 'is_metal'}
        self.calcs_reversed_keys = {'dir_name', 'run_type', 'elements', 'nelements',
                                    'formula_pretty', 'formula_reduced_abc',
                                    'composition_reduced', 'vasp_version', 'formula_anonymous',
                                    'nsites', 'composition_unit_cell', 'completed_at',
                                    'task', 'input', 'output', 'has_vasp_completed'}
        self.analysis_keys = {'delta_volume_percent', 'delta_volume', 'max_force', 'errors',
                              'warnings'}
        self.all_keys = {"root": self.root_keys,
                         "input": self.input_keys,
                         "output": self.output_keys,
                         "calcs_reversed": self.calcs_reversed_keys,
                         "analysis": self.analysis_keys}

    def assimilate(self, path):
        """
        Adapted from matgendb.creator
        Parses vasp runs(vasprun.xml file) and insert the result into the db.

        Args:
            path (str): Path to the directory containing vasprun.xml file

        Returns:
            If in simulate_mode, the entire doc is returned for debugging
            purposes. Else, only the task_id of the inserted doc is returned.
        """
        return self.assimilate_return_task_doc(path)[0]

    def assimilate_return_task_doc(self, path):
        """
        Parses vasp runs(vasprun.xml file) and insert the result into the db.

        Args:
            path (str): Path to the directory containing vasprun.xml file

        Returns:
            If successful, tuple of (task_id, task_doc dict)
            Else, tuple of (False, None)
        """
        try:
            d = self.get_task_doc(path)
            if self.mapi_key is not None and d["state"] == "successful":
                self.calculate_stability(d)
            tid = self._insert_doc(d)
            return tid, d

        except:
            import traceback
            logger.error(traceback.format_exc())
            return False, None

    def get_task_doc(self, path):
        """
        Adapted from matgendb.creator
        Get the entire task doc from the vasprum.xml and the OUTCAR files in the path.
        Also adds some post-processed info.

        Args:
            path (str): Path to the directory containing vasprun.xml and OUTCAR files

        Returns:
            The dictionary to be inserted into the db
        """
        logger.info("Getting task doc for base dir :{}".format(path))
        vasprun_files = self.filter_files(path, file_pattern="vasprun.xml")
        outcar_files = self.filter_files(path, file_pattern="OUTCAR")
        if len(vasprun_files) > 0 and len(outcar_files) > 0:
            d = self.generate_doc(path, vasprun_files, outcar_files)
            self.post_process(path, d)
        else:
            raise ValueError("No VASP files found!")
        self.check_keys(d)
        return d

    def filter_files(self, path, file_pattern="vasprun.xml"):
        """
        Find the files that match the pattern in the given path and
        return them in an ordered dictionary. The searched for files are
        filtered by the run types defined in self.runs. e.g. ["relax1", "relax2"].
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
            fullpath = os.path.abspath(dir_name)
            d = {k: v for k, v in self.additional_fields.items()}
            d["schema"] = {"code": "MatMethods",
                           "version": VaspToDbTaskDrone.__version__}
            d["dir_name"] = fullpath
            d["calcs_reversed"] = [self.process_vasprun(dir_name, taskname, filename)
                                   for taskname, filename in vasprun_files.items()]
            outcar_data = [self.process_outcar(dir_name, filename)
                           for taskname, filename in outcar_files.items()]
            run_stats = {}
            # set run_stats and calcs_reversed.x.output.outcar
            for i, d_calc in enumerate(d["calcs_reversed"]):
                run_stats[d_calc["task"]["name"]] = outcar_data[i].pop("run_stats")
                if d_calc.get("output", {}):
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
            # reverse the calculations data order
            d["calcs_reversed"].reverse()
            d_calc_initial = d["calcs_reversed"][-1]
            d_calc_final = d["calcs_reversed"][0]
            d["chemsys"] = "-".join(sorted(d_calc_final["elements"]))
            comp = Composition(d_calc_final["composition_unit_cell"])
            d["formula_anonymous"] = comp.anonymized_formula
            d["formula_reduced_abc"] = comp.reduced_composition.alphabetical_formula
            for root_key in ["completed_at", "nsites",
                             "composition_unit_cell",
                             "composition_reduced", "formula_pretty",
                             "elements", "nelements"]:
                d[root_key] = d_calc_final[root_key]
            # set other root keys
            self.set_input_data(d_calc_initial, d)
            self.set_output_data(d_calc_final, d)
            self.set_state(d_calc_final, d)
            self.set_analysis(d)
            d["last_updated"] = datetime.datetime.today()
            return d
        except Exception as ex:
            import traceback
            logger.error(traceback.format_exc())
            logger.error("Error in " + os.path.abspath(dir_name) + ".\n" + traceback.format_exc())
            return None

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
        for k, v in {"formula_pretty": "pretty_formula",
                     "composition_reduced": "reduced_cell_formula",
                     "composition_unit_cell": "unit_cell_formula"}.items():
            d[k] = d.pop(v)

        comp = Composition(d["composition_unit_cell"])
        d["formula_anonymous"] = comp.anonymized_formula
        d["formula_reduced_abc"] = comp.reduced_composition.alphabetical_formula
        d["dir_name"] = os.path.abspath(dir_name)
        d["completed_at"] = str(datetime.datetime.fromtimestamp(os.path.getmtime(vasprun_file)))
        d["density"] = vrun.final_structure.density
        # replace 'crystal' with 'structure'
        d["input"]["structure"] = d["input"].pop("crystal")
        d["output"]["structure"] = d["output"].pop("crystal")
        for k, v in {"energy": "final_energy",
                     "energy_per_atom": "final_energy_per_atom"}.items():
            d["output"][k] = d["output"].pop(v)
        if self.parse_dos and self.parse_dos != 'final':
            try:
                d["dos"] = vrun.complete_dos.as_dict()
            except Exception:
                logger.error("No valid dos data exist in {}.\n Skipping dos".format(dir_name))
        if self.bandstructure_mode:
            try:
                bs = vrun.get_band_structure(line_mode=(self.bandstructure_mode == "line"))
                d["bandstructure"] = bs.as_dict()
            except Exception:
                logger.error("No valid band structure data exist in {}.\n Skipping band structure"
                             .format(dir_name))
        else:
            bs = vrun.get_band_structure()
        d["output"]["vbm"] = bs.get_vbm()["energy"]
        d["output"]["cbm"] = bs.get_cbm()["energy"]
        bs_gap = bs.get_band_gap()
        d["output"]["bandgap"] = bs_gap["energy"]
        d["output"]["is_gap_direct"] = bs_gap["direct"]
        d["output"]["is_metal"] = bs.is_metal()
        d["task"] = {"type": taskname, "name": taskname}
        return d

    def process_outcar(self, dir_name, filename):
        """
        Process the outcar file
        """
        outcar_file = os.path.join(dir_name, filename)
        outcar = Outcar(outcar_file)
        return outcar.as_dict()

    def set_input_data(self, d_calc, d):
        """
        set the 'input' key
        """
        # store any overrides to the exchange correlation functional
        xc = d_calc["input"]["incar"].get("GGA")
        if xc:
            xc = xc.upper()
        p = d_calc["input"]["potcar_type"][0].split("_")
        pot_type = p[0]
        functional = "lda" if len(pot_type) == 1 else "_".join(p[1:])
        d["input"] = {"structure": d_calc["input"]["structure"],
                      "is_hubbard": d_calc.pop("is_hubbard"),
                      "hubbards": d_calc.pop("hubbards"),
                      "is_lasph": d_calc["input"]["incar"].get("LASPH", False),
                      "potcar_spec": d_calc["input"].get("potcar_spec"),
                      "xc_override": xc,
                      "pseudo_potential": {"functional": functional.lower(),
                                           "pot_type": pot_type.lower(),
                                           "labels": d_calc["input"]["potcar"]}
                      }

    def set_output_data(self, d_calc, d):
        """
        set the 'output' key
        """
        d["output"] = {
            "structure": d_calc["output"]["structure"],
            "density": d_calc.pop("density"),
            "energy": d_calc["output"]["energy"],
            "energy_per_atom": d_calc["output"]["energy_per_atom"]}
        d["output"].update(self.get_basic_processed_data(d))
        sg = SpacegroupAnalyzer(
            Structure.from_dict(d_calc["output"]["structure"]), 0.1)
        d["output"]["spacegroup"] = {
            "source": "spglib",
            "symbol": sg.get_spacegroup_symbol(),
            "number": sg.get_spacegroup_number(),
            "point_group": sg.get_point_group(),
            "crystal_system": sg.get_crystal_system(),
            "hall": sg.get_hall()}

    def set_state(self, d_calc, d):
        """
        set the 'state' key
        """
        d["state"] = "successful" if d_calc["has_vasp_completed"] else "unsuccessful"

    def set_analysis(self, d, max_force_threshold=0.5, volume_change_threshold=0.2):
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
        if abs(percent_delta_vol) > volume_change_threshold:
            warning_msgs.append("Volume change > {}%"
                                .format(volume_change_threshold * 100))
        max_force = None
        calc = d["calcs_reversed"][0]
        if d["state"] == "successful" and calc["input"]["parameters"].get("NSW", 0) > 0:
            # handle the max force and max force error
            max_force = max([np.linalg.norm(a)
                             for a in calc["output"]
                             ["ionic_steps"][-1]["forces"]])
            if max_force > max_force_threshold:
                error_msgs.append("Final max force exceeds {} eV"
                                  .format(max_force_threshold))
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

    def calculate_stability(self, d):
        """
        Compute the stability of the new entry wrt the existing entries of
        similar composition from the materialsproject database.

        Adapted from pymatgen-db
        """
        m = MPRester(self.mapi_key)
        functional = d["pseudo_potential"]["functional"]
        syms = ["{} {}".format(functional, l) for l in d["pseudo_potential"]["labels"]]
        entry = ComputedEntry(Composition(d["unit_cell_formula"]), d["output"]["final_energy"],
                              parameters={"hubbards": d["hubbards"], "potcar_symbols": syms})
        data = m.get_stability([entry])[0]
        for k in ("e_above_hull", "decomposes_to"):
            d["analysis"][k] = data[k]

    def get_basic_processed_data(self, d):
        """
        return processed data such as vbm, cbm, gap etc.
        """
        calc = d["calcs_reversed"][0]
        gap = calc["output"]["bandgap"]
        cbm = calc["output"]["cbm"]
        vbm = calc["output"]["vbm"]
        is_direct = calc["output"]["is_gap_direct"]
        is_metal = calc["output"]["is_metal"]
        return {
            "bandgap": gap,
            "cbm": cbm,
            "vbm": vbm,
            "is_gap_direct": is_direct,
            "is_metal": is_metal}

    def post_process(self, dir_name, d):
        """
        Simple post-processing for various files other than the vasprun.xml and OUTCAR.
        Looks for files: Transformations.json and custodian.json. Modify this if other
        output files need to be processed.

        Args:
            dir_name:
                The dir_name.
            d:
                Current doc generated.
        """
        logger.info("Post-processing dir:{}".format(dir_name))
        fullpath = os.path.abspath(dir_name)
        # VASP input generated by pymatgen's alchemy has a
        # transformations.json file that keeps track of the origin of a
        # particular structure. This is extremely useful for tracing back a
        # result. If such a file is found, it is inserted into the task doc
        # as d["transformations"]
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
            with zopen(filenames[0], "r") as f:
                d["custodian"] = json.load(f)
        # Convert to full uri path.
        if self.use_full_uri:
            d["dir_name"] = get_uri(dir_name)
        if new_tags:
            d["tags"] = new_tags
        logger.info("Post-processed " + fullpath)

    def check_keys(self, d):
        """
        Sanity check.
        Make sure all the important keys are set
        """
        for k, v in self.all_keys.items():
            if k == "calcs_reversed":
                diff = v.difference(set(d.get(k, d)[0].keys()))
            else:
                diff = v.difference(set(d.get(k, d).keys()))
            if diff:
                logger.warn("The keys {0} in {1} not set".format(diff, k))

    def _insert_doc(self, d):
        if not self.simulate:
            # Perform actual insertion into db. Because db connections cannot
            # be pickled, every insertion needs to create a new connection
            # to the db.
            conn = MongoClient(self.host, self.port)
            db = conn[self.database]
            if self.user:
                db.authenticate(self.user, self.password)
            coll = db[self.collection]
            # Insert dos data into gridfs and then remove it from the dict.
            # DOS data tends to be above the 4Mb limit for mongo docs. A ref
            # to the dos file is in the dos_fs_id.
            result = coll.find_one({"dir_name": d["dir_name"]},
                                   ["dir_name", "task_id"])
            if result is None or self.update_duplicates:
                if self.parse_dos and "calcs_reversed" in d:
                    if "dos" in d["calcs_reversed"][0]:
                        dos = json.dumps(d["calcs_reversed"][0]["dos"],
                                         cls=MontyEncoder)
                        if self.compress_dos:
                            dos = zlib.compress(dos, self.compress_dos)
                            d["calcs_reversed"][0]["dos_compression"] = "zlib"
                        fs = gridfs.GridFS(db, "dos_fs")
                        dosid = fs.put(dos)
                        d["calcs_reversed"][0]["dos_fs_id"] = dosid
                        del d["calcs_reversed"][0]["dos"]
                if self.bandstructure_mode and "calcs_reversed" in d:
                    if "bandstructure" in d["calcs_reversed"][0]:
                        bs = json.dumps(d["calcs_reversed"][0]["bandstructure"], cls=MontyEncoder)
                        if self.compress_bs:
                            bs = zlib.compress(bs, self.compress_bs)
                            d["calcs_reversed"][0]["bandstructure_compression"] = "zlib"
                        fs = gridfs.GridFS(db, "bandstructure_fs")
                        bsid = fs.put(bs)
                        d["calcs_reversed"][0]["bandstructure_fs_id"] = bsid
                        del d["calcs_reversed"][0]["bandstructure"]
                d["last_updated"] = datetime.datetime.today()
                if result is None:
                    if ("task_id" not in d) or (not d["task_id"]):
                        d["task_id"] = db.counter.find_and_modify(
                            query={"_id": "taskid"},
                            update={"$inc": {"c": 1}}
                        )["c"]
                    logger.info("Inserting {} with taskid = {}"
                                .format(d["dir_name"], d["task_id"]))
                elif self.update_duplicates:
                    d["task_id"] = result["task_id"]
                    logger.info("Updating {} with taskid = {}"
                                .format(d["dir_name"], d["task_id"]))

                coll.update({"dir_name": d["dir_name"]}, {"$set": d},
                            upsert=True)
                return d["task_id"]
            else:
                logger.info("Skipping duplicate {}".format(d["dir_name"]))
        else:
            d["task_id"] = 0
            logger.info("Simulated insert into database for {} with task_id {}"
                        .format(d["dir_name"], d["task_id"]))
            return d

    def get_valid_paths(self, path):
        """
        Required by the AbstractDrone.
        Update this and use it to further filter the files to be assimilated.
        """
        pass

    def as_dict(self):
        init_args = {"host": self.host, "port": self.port,
                     "database": self.database,
                     "user": self.user,
                     "password": self.password,
                     "collection": self.collection,
                     "parse_dos": self.parse_dos,
                     "compress_dos": self.compress_dos,
                     "bandstructure_mode": self.bandstructure_mode,
                     "compress_bs": self.compress_bs,
                     "simulate_mode": self.simulate,
                     "additional_fields": self.additional_fields,
                     "update_duplicates": self.update_duplicates,
                     "mapi_key": self.mapi_key,
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

    @classmethod
    def from_db_doc(cls, dbdoc=None, additional_fields=None, options=None):
        """
        Set the drone from the database connection settings.
        """
        additional_fields = additional_fields if additional_fields else {}
        options = options if options else {}

        if dbdoc:
            return VaspToDbTaskDrone(host=dbdoc["host"], port=dbdoc["port"],
                                     database=dbdoc["database"],
                                     user=dbdoc.get("admin_user"),
                                     password=dbdoc.get("admin_password"),
                                     collection=dbdoc["collection"],
                                     additional_fields=additional_fields,
                                     parse_dos=options.get("parse_dos", False),
                                     compress_dos=options.get("compress_dos", False),
                                     bandstructure_mode=options.get("bandstructure_mod", False),
                                     compress_bs=options.get("compress_bs", False),
                                     update_duplicates=options.get("update_duplicates", True),
                                     mapi_key=options.get("mapi_key", None),
                                     use_full_uri=options.get("use_full_uri", True),
                                     runs=options.get("runs", ["relax1", "relax2"]))
        else:
            return VaspToDbTaskDrone(simulate_mode=True)
