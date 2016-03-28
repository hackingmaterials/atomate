# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
This module defines the drones
"""

import os
import string
import datetime
from fnmatch import fnmatch
from collections import OrderedDict

from monty.os.path import zpath

import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Vasprun

from matgendb.creator import VaspToDbTaskDrone


__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'
__date__ = 'Mar 27, 2016'


class MMVaspToDbTaskDrone(VaspToDbTaskDrone):
    __version__ = 0.1

    def __init__(self, host="127.0.0.1", port=27017, database="vasp",
                 user=None, password=None, collection="tasks",
                 parse_dos=False, compress_dos=False, simulate_mode=False,
                 additional_fields=None, update_duplicates=True,
                 mapi_key=None, use_full_uri=True, runs=None):
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

        Parses vasp runs and insert the result into the db.

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

        Get the entire task doc for a path, including any post-processing.
        """
        print("Getting task doc for base dir :{}".format(path))
        files = os.listdir(path)
        vasprun_files = OrderedDict()
        if "STOPCAR" in files:
            print(path + " contains stopped run")
        for f in files:  # get any vasprun from the folder
            if fnmatch(f, "vasprun.xml*") and \
                            f not in vasprun_files.values():
                vasprun_files['standard'] = f
        d = {}
        if len(vasprun_files) > 0:
            d = self.generate_doc(path, vasprun_files)
            self.post_process(path, d)
        else:
            raise ValueError("No VASP files found!")
        return d

    def generate_doc(self, dir_name, vasprun_files):
        """
        Adapted from matgendb.creator.generate_doc
        """
        try:
            fullpath = os.path.abspath(dir_name)
            d = {k: v for k, v in self.additional_fields.items()}
            d["name"] = "MatMethods"
            d["dir_name"] = fullpath
            d["schema_version"] = MMVaspToDbTaskDrone.__version__
            d["calculations"] = [
                self.process_vasprun(dir_name, taskname, filename)
                for taskname, filename in vasprun_files.items()]
            d1 = d["calculations"][0]
            d2 = d["calculations"][-1]
            d["chemsys"] = "-".join(sorted(d2["elements"]))
            vals = sorted(d2["reduced_cell_formula"].values())
            d["anonymous_formula"] = {string.ascii_uppercase[i]: float(vals[i])
                                      for i in range(len(vals))}
            for root_key in ["completed_at", "nsites",
                             "unit_cell_formula",
                             "reduced_cell_formula", "pretty_formula",
                             "elements", "nelements", "run_type"]:
                d[root_key] = d2[root_key]
            self.set_input_data(d1, d2, d)
            self.set_output_data(d1, d2, d)
            self.set_state(vasprun_files, d1, d2, d)
            self.set_analysis(d)
            d["last_updated"] = datetime.datetime.today()
            return d
        except Exception as ex:
            import traceback
            print(traceback.format_exc())
            print("Error in " + os.path.abspath(dir_name) +
                  ".\n" + traceback.format_exc())
            return None

    def process_vasprun(self, dir_name, taskname, filename):
        """
        Adapted from matgendb.creator

        Process a vasprun.xml file.
        """
        vasprun_file = os.path.join(dir_name, filename)
        r = Vasprun(vasprun_file)
        d = r.as_dict()
        d["dir_name"] = os.path.abspath(dir_name)
        d["completed_at"] = \
            str(datetime.datetime.fromtimestamp(os.path.getmtime(
                vasprun_file)))
        d["density"] = r.final_structure.density
        # replace 'crystal' with 'structure'
        d["input"]["structure"] = d["input"].pop("crystal")
        d["output"]["structure"] = d["output"].pop("crystal")
        if self.parse_dos and (self.parse_dos != 'final' \
                                       or taskname == self.runs[-1]):
            try:
                d["dos"] = r.complete_dos.as_dict()
            except Exception:
                print("No valid dos data exist in {}.\n Skipping dos"
                      .format(dir_name))
        d["task"] = {"type": taskname, "name": taskname}
        return d

    def set_input_data(self, d1, d2, d):
        """
        set the 'input' key
        """
        # store any overrides to the exchange correlation functional
        xc = d2["input"]["incar"].get("GGA")
        if xc:
            xc = xc.upper()
        p = d2["input"]["potcar_type"][0].split("_")
        pot_type = p[0]
        functional = "lda" if len(pot_type) == 1 else "_".join(p[1:])
        d["input"] = {"structure": d1["input"]["structure"],
                      "is_hubbard": d2["is_hubbard"],
                      "hubbards": d2["hubbards"],
                      "is_lasph": d2["input"]["incar"].get("LASPH", False),
                      "potcar_spec": d1["input"].get("potcar_spec"),
                      "xc_override": xc,
                      "pseudo_potential": {"functional": functional.lower(),
                                           "pot_type": pot_type.lower(),
                                           "labels": d2["input"]["potcar"]}
                      }

    def set_output_data(self, d1, d2, d):
        """
        set the 'output' key
        """
        d["output"] = {
            "structure": d2["output"]["structure"],
            "density": d2["density"],
            "final_energy": d2["output"]["final_energy"],
            "final_energy_per_atom": d2["output"]["final_energy_per_atom"]
        }
        d["output"].update(self.get_basic_processed_data(d))
        sg = SpacegroupAnalyzer(Structure.from_dict(d2["output"]["structure"]),
                                0.1)
        d["output"]["spacegroup"] = {
            "source": "spglib",
            "symbol": sg.get_spacegroup_symbol(),
            "number": sg.get_spacegroup_number(),
            "point_group": sg.get_point_group(),
            "crystal_system": sg.get_crystal_system(),
            "hall": sg.get_hall()}

    def set_state(self, vasprun_files, d1, d2, d):
        """
        set the 'state' key
        """
        if len(d["calculations"]) == len(self.runs) or \
                        list(vasprun_files.keys())[0] != "relax1":
            d["state"] = "successful" if d2["has_vasp_completed"] \
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
                        d["calculations"][0]["input"]["parameters"].get("NSW",
                                                                        0) > 0:
            # handle the max force and max force error
            max_force = max([np.linalg.norm(a)
                             for a in d["calculations"][-1]["output"]
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
        calc = d["calculations"][-1]
        gap = calc["output"]["bandgap"]
        cbm = calc["output"]["cbm"]
        vbm = calc["output"]["vbm"]
        is_direct = calc["output"]["is_gap_direct"]
        return {
            "bandgap": gap,
            "cbm": cbm,
            "vbm": vbm,
            "is_gap_direct": is_direct}


def is_valid_vasp_dir(mydir):
    """
    Copied from matgendb.creator

    note: that the OUTCAR and POSCAR are known to be empty in some situations
    """
    files = ["OUTCAR", "POSCAR", "INCAR", "KPOINTS"]
    for f in files:
        m_file = os.path.join(mydir, f)
        if not os.path.exists(zpath(m_file)) or not (
                        os.stat(m_file).st_size > 0 or os.stat(
                        m_file + '.gz').st_size > 0):
            return False
    return True
