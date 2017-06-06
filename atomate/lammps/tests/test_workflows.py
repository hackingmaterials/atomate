# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import json
import os
import shutil
import unittest
import filecmp
from pymongo import MongoClient

from pymatgen.io.lammps.data import LammpsForceFieldData

from fireworks import LaunchPad, FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.lammps.workflows.core import wf_from_input_template
from atomate.utils.testing import AtomateTest

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "common", "test_files")

DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test
LAMMPS_CMD = None  # "lmp_serial"


class TestLammpsWorkflows(AtomateTest):

    def setUp(self):
        super(TestLammpsWorkflows, self).setUp()
        self.data_file = os.path.join(module_dir, "test_files/lammps_data.dat")
        self.input_template = os.path.join(module_dir, "test_files/NPT.json")
        self.reference_files_path = os.path.abspath(os.path.join(module_dir, "reference_files"))

    def test_lammps_wflow(self):
        lammps_data = LammpsForceFieldData.from_file(self.data_file)
        log_filename = "peo.log"
        dump_filename = "peo.dump"
        dcd_traj_filename = "peo.dcd"
        timestep = 1  # in fmsec for 'real' units
        n_timesteps = 1000
        dump_freq = 50
        T = [300, 300, 100.0]  # start, end, damp
        P = [0, 0, 100.0]

        # override default settings read from the json file with these
        user_settings = {"log": log_filename,
                         "timestep": timestep,
                         "run": n_timesteps,
                         "pair_style": "buck/coul/cut 15.0",
                         "pair_coeff": ["1 1 2649.6 0.2674 27.22",
                                        "1 2 4320.0 0.2928 137.6",
                                        "1 3 14176.0 0.2563 104.0",
                                        "2 2 14976.0 0.3236 637.6",
                                        "2 3 33702.4 0.2796 503.0",
                                        "3 3 75844.8 0.2461 396.9"],
                         "thermo_style": "custom step time temp press pe ke etotal enthalpy fmax fnorm",
                         "fix": "NPT all npt temp {T[0]} {T[1]} {T[2]} iso {P[0]} {P[1]} {P[2]}".format(T=T, P=P),
                         "dump": [
                             "peodump all custom {} {} id type x y z ix iy iz mol".format(dump_freq, dump_filename),
                             "traj all dcd {} {}".format(dump_freq, dcd_traj_filename)]}

        if not LAMMPS_CMD:
            # fake run
            lammps_bin = "cp  {}/peo.* .".format(self.reference_files_path)
            dry_run = True
        else:
            lammps_bin = LAMMPS_CMD
            dry_run = False
        wf = wf_from_input_template(self.input_template, lammps_data, "npt.data", user_settings,
                                    is_forcefield=True, input_filename="lammps.inp",
                                    lammps_bin=lammps_bin, db_file=">>db_file<<", dry_run=dry_run)
        self.lp.add_wf(wf)
        # run
        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))
        d = self.get_task_collection().find_one()
        self._check_run(d)

    def _check_run(self, d):
        if not LAMMPS_CMD:
            self.assertTrue(filecmp.cmp(os.path.join(d["dir_name"], "peo.log"),
                                        "{}/peo.log".format(self.reference_files_path)))
            self.assertTrue(filecmp.cmp(os.path.join(d["dir_name"], "peo.dump"),
                                        "{}/peo.dump".format(self.reference_files_path)))
            self.assertTrue(filecmp.cmp(os.path.join(d["dir_name"], "peo.dcd"),
                                        "{}/peo.dcd".format(self.reference_files_path)))


if __name__ == "__main__":
    unittest.main()
