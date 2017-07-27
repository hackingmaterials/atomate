# coding: utf-8
from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
import filecmp

from fireworks import FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.lammps.firetasks.run_calc import RunLammpsFake
from atomate.lammps.workflows.core import get_wf_from_input_template
from atomate.utils.testing import AtomateTest

__author__ = 'Kiran Mathew, Brandon Wood'
__email__ = 'kmathew@lbl.gov, b.wood@berkeley.edu'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "common", "test_files")

LAMMPS_CMD = None  # "mpirun -n 4 lmp_mpi"


def use_fake_lammps(original_wf, ref_dir):
    for idx_fw, fw in enumerate(original_wf.fws):
        for idx_t, t in enumerate(fw.tasks):
            if "RunLammps" in str(t):
                original_wf.fws[idx_fw].tasks[idx_t] = RunLammpsFake(ref_dir=ref_dir)
    return original_wf


class TestLammpsWorkflows(AtomateTest):

    def setUp(self):
        super(TestLammpsWorkflows, self).setUp()
        self.data_file = os.path.join(module_dir, "test_files/peo.data")
        self.input_file_template = os.path.join(module_dir, "test_files/peo.in.template")
        self.user_settings = {"pair_style": "buck/coul/long 15.0",
                              "kspace_style": "ewald 1.0e-5",
                              "thermo_style_1": "custom step temp press vol density pe ke etotal enthalpy fmax fnorm",
                              "thermo_style_2": "100",
                              "fix_1": "CHbonds CH shake 0.0001 20 0 b 2",
                              "fix_2": "NPT all npt temp 350 350 100.0 iso 1.0 1.0 1500.0"}
        self.input_filename = "lammps.in"
        self.db_file = os.path.join(db_dir, "db.json")
        self.reference_files_path = os.path.abspath(os.path.join(module_dir, "test_files"))

    def test_lammps_wflow(self):

        wf = get_wf_from_input_template(self.input_file_template, self.user_settings,
                                        lammps_data=self.data_file, is_forcefield=True,
                                        input_filename=self.input_filename,
                                        data_filename="lammps.data", db_file=self.db_file,
                                        name="peo_test")

        if not LAMMPS_CMD:
            wf = use_fake_lammps(wf, self.reference_files_path)

        self.lp.add_wf(wf)

        # run
        rapidfire(self.lp, fworker=FWorker(env={"db_file": self.db_file}))

        d = self.get_task_collection().find_one()

        self._check_run(d)

    def _check_run(self, d):
        if not LAMMPS_CMD:
            path = d["dir_name"].split(":")[-1]
            self.assertTrue(filecmp.cmp(os.path.join(path, "peo.in"),
                                        "{}/peo.in".format(self.reference_files_path)))
            self.assertTrue(filecmp.cmp(os.path.join(path, "lammps.log"),
                                        "{}/lammps.log".format(self.reference_files_path)))
            self.assertTrue(filecmp.cmp(os.path.join(path, "peo.dump"),
                                        "{}/peo.dump".format(self.reference_files_path)))
            self.assertTrue(filecmp.cmp(os.path.join(path, "peo.dcd"),
                                        "{}/peo.dcd".format(self.reference_files_path)))

if __name__ == "__main__":
    unittest.main()
