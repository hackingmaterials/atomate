# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
This module defines the drones
"""

import os
import sys
import re
import string
import datetime
import zlib
from fnmatch import fnmatch
from collections import OrderedDict
import json
import glob
import logging

from monty.io import zopen
from monty.json import MontyEncoder

import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Vasprun, Outcar

from pymongo import MongoClient
import gridfs

from matgendb.creator import VaspToDbTaskDrone, get_uri

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'
__date__ = 'Mar 27, 2016'

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
sh = logging.StreamHandler(stream=sys.stdout)
sh.setFormatter(formatter)
logger.addHandler(sh)


class MMVaspToDbTaskDrone(VaspToDbTaskDrone):
    """
    VaspToDbTaskDrone with updated schema.
    Also removed the processing of aflow style runs.
    Please refer to matgendb.creator.VaspToDbTaskDrone documentation
    """

    __version__ = 0.1

    def __init__(self, host="127.0.0.1", port=27017,
                 database="vasp", collection="tasks",
                 user=None, password=None, taskname="standard",
                 parse_dos=False, compress_dos=False, simulate_mode=False,
                 additional_fields=None, update_duplicates=True,
                 mapi_key=None, use_full_uri=True, runs=None):
        self.taskname = taskname
        self.root_keys = {"name", "dir_name", "schema_version", "chemsys",
                          "anonymous_formula", "calculations", "completed_at",
                          "nsites", "unit_cell_formula",
                          "reduced_cell_formula", "pretty_formula",
                          "elements", "nelements", "run_type",
                          "input", "output", "state", "analysis"}
        self.input_keys = {'is_lasph', 'is_hubbard', 'xc_override',
                           'potcar_spec', 'hubbards', 'structure',
                           'pseudo_potential'}
        self.output_keys = {'is_gap_direct', 'density', 'bandgap',
                            'final_energy_per_atom', 'vbm', 'cbm',
                            'spacegroup', 'final_energy', 'structure'}
        self.calculations_keys = {'dir_name', 'run_type', 'elements',
                                  'hubbards', 'nelements', 'pretty_formula',
                                  'reduced_cell_formula', 'vasp_version',
                                  'nsites', 'unit_cell_formula',
                                  'completed_at', 'output',
                                  'task', 'is_hubbard', 'input', 'task',
                                  'has_vasp_completed'}
        self.analysis_keys = {'delta_volume_percent', 'delta_volume',
                              'max_force', 'errors', 'warnings'}
        self.all_keys = {"root": self.root_keys, "input": self.input_keys,
                         "output": self.output_keys,
                         "calculations": self.calculations_keys,
                         "analysis": self.analysis_keys}
        super(MMVaspToDbTaskDrone, self).__init__(host=host, port=port,
                                                  database=database,
                                                  user=user, password=password,
                                                  collection=collection,
                                                  parse_dos=parse_dos,
                                                  compress_dos=compress_dos,
                                                  simulate_mode=simulate_mode,
                                                  additional_fields=additional_fields,
                                                  update_duplicates=update_duplicates,
                                                  mapi_key=mapi_key,
                                                  use_full_uri=use_full_uri,
                                                  runs=runs
                                                  )

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
        try:
            d = self.get_task_doc(path)
            if self.mapi_key is not None and d["state"] == "successful":
                self.calculate_stability(d)
            tid = self._insert_doc(d)
            return tid
        except Exception as ex:
            import traceback
            logger.error(traceback.format_exc())
            return False

    def get_task_doc(self, path):
        """
        Adapted from matgendb.creator
        Get the entire task doc from the vasprum.xml file in the path.
        Processes only one xml file from the given path.
        Also adds some post-processed info.

        Args:
            path (str): Path to the directory containing vasprun.xml file

        Returns:
            The dictionary to be inserted into the db
        """
        logger.info("Getting task doc for base dir :{}".format(path))
        files = os.listdir(path)
        vasprun_file = OrderedDict()
        if "STOPCAR" in files:
            logger.warn(path + " contains stopped run")
        for f in files:
            if fnmatch(f, "vasprun.xml*"):
                vasprun_file['standard'] = f
        d = {}
        if vasprun_file:
            d = self.generate_doc(path, vasprun_file)
            self.post_process(path, d)
        else:
            raise ValueError("No VASP files found!")
        self.check_keys(d)
        return d

    def generate_doc(self, dir_name, vasprun_file):
        """
        Adapted from matgendb.creator.generate_doc
        """
        try:
            fullpath = os.path.abspath(dir_name)
            d = {k: v for k, v in self.additional_fields.items()}
            d["name"] = "MatMethods"
            d["dir_name"] = fullpath
            d["schema_version"] = MMVaspToDbTaskDrone.__version__
            filename = vasprun_file[self.taskname]
            d["calculations"] = self.process_vasprun(dir_name, filename)
            d_calc = d["calculations"]
            d["chemsys"] = "-".join(sorted(d_calc["elements"]))
            vals = sorted(d_calc["reduced_cell_formula"].values())
            d["anonymous_formula"] = {string.ascii_uppercase[i]: float(vals[i])
                                      for i in range(len(vals))}
            for root_key in ["completed_at", "nsites",
                             "unit_cell_formula",
                             "reduced_cell_formula", "pretty_formula",
                             "elements", "nelements", "run_type"]:
                d[root_key] = d_calc[root_key]
            self.set_input_data(d_calc, d)
            self.set_output_data(d_calc, d)
            self.set_state(vasprun_file, d_calc, d)
            self.set_analysis(d)
            d["last_updated"] = datetime.datetime.today()
            return d
        except Exception as ex:
            import traceback
            logger.error(traceback.format_exc())
            logger.error("Error in " + os.path.abspath(dir_name) +
                         ".\n" + traceback.format_exc())
            return None

    def process_vasprun(self, dir_name, filename):
        """
        Adapted from matgendb.creator

        Process a vasprun.xml file.
        """
        vasprun_file = os.path.join(dir_name, filename)
        vrun = Vasprun(vasprun_file)
        d = vrun.as_dict()
        d["dir_name"] = os.path.abspath(dir_name)
        d["completed_at"] = \
            str(datetime.datetime.fromtimestamp(os.path.getmtime(
                vasprun_file)))
        d["density"] = vrun.final_structure.density
        # replace 'crystal' with 'structure'
        d["input"]["structure"] = d["input"].pop("crystal")
        d["output"]["structure"] = d["output"].pop("crystal")
        if self.parse_dos and self.parse_dos != 'final':
            try:
                d["dos"] = vrun.complete_dos.as_dict()
            except Exception:
                logger.error("No valid dos data exist in {}.\n Skipping dos"
                             .format(dir_name))
        d["task"] = {"type": self.taskname, "name": self.taskname}
        return d

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
                      "is_hubbard": d_calc["is_hubbard"],
                      "hubbards": d_calc["hubbards"],
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
            "final_energy": d_calc["output"]["final_energy"],
            "final_energy_per_atom": d_calc["output"]["final_energy_per_atom"]}
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

    def set_state(self, vasprun_file, d_calc, d):
        """
        set the 'state' key
        """
        if list(vasprun_file.keys())[0] == "standard":
            d["state"] = "successful" if d_calc["has_vasp_completed"] \
                else "unsuccessful"
        else:
            d["state"] = "stopped"

    def set_analysis(self, d, max_force_threshold=0.5,
                     volume_change_threshold=0.2):
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
        if d["state"] == "successful" and \
                        d["calculations"]["input"]["parameters"].get("NSW",
                                                                     0) > 0:
            # handle the max force and max force error
            max_force = max([np.linalg.norm(a)
                             for a in d["calculations"]["output"]
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

    def get_basic_processed_data(self, d):
        """
        return processed data such as vbm, cbm, gap
        """
        calc = d["calculations"]
        gap = calc["output"]["bandgap"]
        cbm = calc["output"]["cbm"]
        vbm = calc["output"]["vbm"]
        is_direct = calc["output"]["is_gap_direct"]
        return {
            "bandgap": gap,
            "cbm": cbm,
            "vbm": vbm,
            "is_gap_direct": is_direct}

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
                if self.parse_dos and "calculations" in d:
                    if "dos" in d["calculations"]:
                        dos = json.dumps(d["calculations"]["dos"],
                                         cls=MontyEncoder)
                        if self.compress_dos:
                            dos = zlib.compress(dos, self.compress_dos)
                            d["calculations"]["dos_compression"] = "zlib"
                        fs = gridfs.GridFS(db, "dos_fs")
                        dosid = fs.put(dos)
                        d["calculations"]["dos_fs_id"] = dosid
                        del d["calculations"]["dos"]
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

    def post_process(self, dir_name, d):
        """
        Simple post-processing for various files other than the vasprun.xml.
        Called by generate_task_doc. Modify this if your runs have other
        kinds of processing requirements.

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
                    m = re.match("(\d+)-ICSD",
                                 transformations["history"][0]["source"])
                    if m:
                        d["icsd_id"] = int(m.group(1))
                except Exception as ex:
                    logger.warning("Cannot parse ICSD from transformations "
                                   "file.")
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
        # Parse OUTCAR for additional information and run stats that are
        # generally not in vasprun.xml.
        try:
            run_stats = {}
            for filename in glob.glob(os.path.join(fullpath, "OUTCAR*")):
                outcar = Outcar(filename)
                d["calculations"]["output"]["outcar"] = outcar.as_dict()
                run_stats[self.taskname] = outcar.run_stats
        except:
            logger.error("Bad OUTCAR for {}.".format(fullpath))
        try:
            overall_run_stats = {}
            for key in ["Total CPU time used (sec)", "User time (sec)",
                        "System time (sec)", "Elapsed time (sec)"]:
                overall_run_stats[key] = sum([v[key]
                                              for v in run_stats.values()])
            run_stats["overall"] = overall_run_stats
        except:
            logger.error("Bad run stats for {}.".format(fullpath))
        d["run_stats"] = run_stats
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
            diff = v.difference(set(d.get(k, d).keys()))
            if diff:
                logger.warn("The keys {0} in {1} not set".format(diff, k))

    @classmethod
    def from_dict(cls, d):
        return cls(**d["init_args"])

    def as_dict(self):
        init_args = {"host": self.host, "port": self.port,
                     "database": self.database, "user": self.user,
                     "password": self.password,
                     "collection": self.collection,
                     "taskname": self.taskname,
                     "parse_dos": self.parse_dos,
                     "simulate_mode": self.simulate,
                     "additional_fields": self.additional_fields,
                     "update_duplicates": self.update_duplicates}
        output = {"name": self.__class__.__name__,
                  "init_args": init_args, "version":
                      self.__class__.__version__}
        return output
