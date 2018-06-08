# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import json
import os
import shutil
import unittest

import numpy as np

from pymongo import MongoClient

from fireworks import LaunchPad, FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.utils.testing import AtomateTest
from atomate.vasp.powerups import use_fake_vasp
from atomate.vasp.workflows.base.ferroelectric import get_wf_ferroelectric
from atomate.utils.utils import get_a_unique_id
from atomate.vasp.firetasks.parse_outputs import PolarizationToDb

from pymatgen import SETTINGS

__author__ = 'Tess Smidt'
__email__ = 'blondegeek@gmail.com'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
ref_dir = os.path.join(module_dir, "..", "..", "test_files")

from pymatgen.core.structure import Structure

DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...

@unittest.skip("This workflow may require a set of POTCARs and corresponding .pmgrc.yaml file.")
class TestFerroelectricWorkflow(AtomateTest):

    def setUp(self):
        super(TestFerroelectricWorkflow, self).setUp()

        self.bto_polar = Structure.from_file(ref_dir+"/ferroelectric_wf/"+"BTO_polar_POSCAR")
        self.bto_nonpolar = Structure.from_file(ref_dir+"/ferroelectric_wf/"+"BTO_nonpolar_POSCAR")

        self.wfid = "wfid_" + get_a_unique_id()

        self.ferroelectric_config = {'vasp_cmd': '>>vasp_cmd<<',
                                    'db_file': '>>db_file<<',
                                    'nimages': 2,
                                    'relax' : True,
                                    'wfid': self.wfid,
                                    'add_analysis_task': False}

        self.wf = get_wf_ferroelectric(self.bto_polar, self.bto_nonpolar, **self.ferroelectric_config)

        if not os.path.exists(os.path.join(ref_dir, "ferroelectric_wf", "nonpolar_static")):
            self.untarTestFiles()

    def untarTestFiles(self):
        import tarfile
        tar_filename = os.path.abspath(os.path.join(ref_dir, "ferroelectric_wf",
                                    "test_ferroelectric_workflow.gz.tar"))
        print(tar_filename)
        t = tarfile.open(tar_filename)
        t.extractall(os.path.abspath(os.path.join(ref_dir, "ferroelectric_wf")))

    def _simulate_vasprun(self, wf):
        reference_dir = os.path.abspath(os.path.join(ref_dir, "ferroelectric_wf"))
        bto_ref_dirs = {"_polar_relaxation": os.path.join(reference_dir, "polar_relaxation"),
                        "_polar_static": os.path.join(reference_dir, "polar_static"),
                        "_polar_polarization": os.path.join(reference_dir, "polar_polarization"),
                        "_nonpolar_relaxation": os.path.join(reference_dir, "nonpolar_relaxation"),
                        "_nonpolar_static": os.path.join(reference_dir, "nonpolar_static"),
                        "_nonpolar_polarization": os.path.join(reference_dir, "nonpolar_polarization"),
                        "_interpolation_1_static": os.path.join(reference_dir, "interpolation_1_static"),
                        "_interpolation_1_polarization": os.path.join(reference_dir, "interpolation_1_polarization")}
        # Add test with for analysis?
        # Which params_to_check?
        return use_fake_vasp(wf, bto_ref_dirs, params_to_check=["ENCUT", "LWAVE"])

    def _check_run(self, d, mode):

        # Check polar and nonpolar relaxations
        if mode == '_polar_relaxation':
            self.assertAlmostEqual(d["calcs_reversed"][0]["output"]["structure"]["lattice"]["c"], 4.2157, 2)

        if mode == '_nonpolar_relaxation':
            self.assertAlmostEqual(d["calcs_reversed"][0]["output"]["structure"]["lattice"]["c"], 4.0350, 2)

        # Check interpolated structure
        if mode == '_interpolation_1_polarization':
            self.assertAlmostEqual(d["calcs_reversed"][0]["output"]["structure"]["lattice"]["c"], 4.1345, 2)

        # Check that Outcar has needed keys for polarization analysis.
        if '_polarization' in mode and 'processing' not in mode:
            # Check that Outcar has p_ion, p_elec, zval_dict
            assert d["calcs_reversed"][0]["output"]["outcar"].get("p_ion", None) is not None
            assert d["calcs_reversed"][0]["output"]["outcar"].get("p_elec", None) is not None
            assert d["calcs_reversed"][0]["output"]["outcar"].get("zval_dict", None) is not None

        # Check analysis step.
        if mode == "_polarization_post_processing":
            self.assertAlmostEqual(d['polarization_change_norm'], 46.288752795325244)

    def test_wf(self):

        self.wf = self._simulate_vasprun(self.wf)

        # 2*relax + 3*polarization = 5
        self.assertEqual(len(self.wf.fws), 5)

        # check VASP parameters on polarization calculation for interpolated structures
        interpolated_polarization_vis = [fw.tasks[7]['incar_update']['lcalcpol']
            for fw in self.wf.fws if "polarization" in fw.name and "interpolation" in fw.name]

        assert all(interpolated_polarization_vis)

        self.lp.add_wf(self.wf)
        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))

        # Check polar relaxation
        d = self.get_task_collection().find_one({"task_label": "_polar_relaxation"})
        self._check_run(d, "_polar_relaxation")

        # Check nonpolar relaxation
        d = self.get_task_collection().find_one({"task_label": "_nonpolar_relaxation"})
        self._check_run(d, "_nonpolar_relaxation")

        # Check polarization calculations
        D = self.get_task_collection().find({"task_label": {"$regex": ".*polarization"}})
        for d in D:
            self._check_run(d, d["task_label"])


if __name__ == "__main__":
    unittest.main()
