# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

import numpy as np

from pymatgen import Structure
from pymatgen.io.feff.inputs import Tags

from fireworks.core.fworker import FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.feff.workflows.core import get_wf_xas
from atomate.utils.testing import AtomateTest

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")

DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test
FEFF_CMD = None  # "feff"


class TestXASWorkflow(AtomateTest):

    def setUp(self):
        super(TestXASWorkflow, self).setUp()
        # CoO
        self.structure = Structure.from_file(os.path.join(module_dir, "..", "..", "test_files", "Co2O2.cif"))
        #PymatgenTest.get_mp_structure("mp-715460")
        self.user_tag_settings = {"RPATH": -1,
                                 "SCF": "7 0 30 0.2 3",
                                 "FMS": "9 0",
                                 "LDOS": "-30.0 30.0 0.1",
                                 "RECIPROCAL": "",
                                 "EDGE": "K"}
        # 3rd site
        self.absorbing_atom = 2
        self.edge = "K"
        self.nkpts = 1000

    def test_xas_wflow_abatom_by_idx(self):
        if not FEFF_CMD:
            # fake run
            xmu_file_path = os.path.abspath(os.path.join(module_dir, "../../test_files/xmu.dat"))
            feff_bin = "cp {} .".format(xmu_file_path)
        else:
            feff_bin = FEFF_CMD

        wf = get_wf_xas(self.absorbing_atom, self.structure, feff_input_set="XANES", edge="K",
                        feff_cmd=feff_bin, db_file=">>db_file<<", use_primitive=False,
                        user_tag_settings=self.user_tag_settings)
        self.assertEqual(len(wf.as_dict()["fws"]), 1)

        self.lp.add_wf(wf)
        # run
        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))

        d = self.get_task_collection().find_one({"spectrum_type": "XANES"})
        self._check_run(d)

    def test_xas_wflow_abatom_by_symbol(self):
        wf_prim = get_wf_xas("O", self.structure, feff_input_set="XANES", edge="K",
                             use_primitive=True, user_tag_settings=self.user_tag_settings)
        wf = get_wf_xas("O", self.structure, feff_input_set="XANES", edge="K",
                        use_primitive=False, user_tag_settings=self.user_tag_settings)
        self.assertEqual(len(wf_prim.as_dict()["fws"]), 1)
        self.assertEqual(len(wf.as_dict()["fws"]), 2)

    def test_xanes_vs_exafs(self):
        wf_xanes = get_wf_xas(self.absorbing_atom, self.structure, feff_input_set="XANES", edge="K",
                              user_tag_settings=self.user_tag_settings)
        wf_exafs = get_wf_xas(self.absorbing_atom, self.structure, feff_input_set="EXAFS", edge="K")

        self.assertEqual(wf_xanes.as_dict()["fws"][0]["spec"]['_tasks'][0]['feff_input_set']['@class'],
                         'MPXANESSet')
        self.assertEqual(wf_exafs.as_dict()["fws"][0]["spec"]['_tasks'][0]['feff_input_set']['@class'],
                         'MPEXAFSSet')

    def _check_run(self, d):
        run_dir = d["dir_name"]
        self.assertEqual(d["edge"], self.edge)
        self.assertEqual(d["absorbing_atom"], self.absorbing_atom)
        tags = Tags.from_file(os.path.join(run_dir, "feff.inp"))
        self.assertEqual(d["input_parameters"], tags.as_dict())
        xmu = np.loadtxt((os.path.join(run_dir, "xmu.dat")))
        self.assertEqual(d["spectrum"], xmu.tolist())


if __name__ == "__main__":
    unittest.main()
