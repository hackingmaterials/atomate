# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

from pymatgen import Structure
from pymatgen.io.feff.inputs import Tags

from fireworks.core.fworker import FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.feff.workflows.core import get_wf_eels
from atomate.utils.testing import AtomateTest

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test


class TestEELSWorkflow(AtomateTest):

    def setUp(self):
        super(TestEELSWorkflow, self).setUp()
        self.structure = Structure.from_file(os.path.join(module_dir, "..", "..", "test_files",
                                                          "Co2O2.cif"))
        self.user_tag_settings = {"RPATH": -1,
                                 "SCF": "7 0 30 0.2 3",
                                 "FMS": "9 0",
                                 "LDOS": "-30.0 30.0 0.1",
                                 "RECIPROCAL":"",
                                 "EDGE": "L1",
                                 "COREHOLE": "RPA"}
        # 3rd site
        self.absorbing_atom = 2
        self.edge = "L1"
        self.nkpts = 1000

    def test_eels_wflow_abatom_by_idx(self):
        # for the sake of test just copy xmu to eels
        xmu_file_path = os.path.abspath(os.path.join(module_dir, "../../test_files/xmu.dat"))
        feff_bin = "cp {} eels.dat".format(xmu_file_path)
        wf = get_wf_eels(self.absorbing_atom, self.structure, feff_input_set="ELNES",
                        edge="L1", user_tag_settings=self.user_tag_settings, use_primitive=False,
                         feff_cmd=feff_bin, db_file=">>db_file<<")
        self.assertEqual(len(wf.as_dict()["fws"]), 1)

        self.lp.add_wf(wf)
        # run
        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))

        d = self.get_task_collection().find_one({"spectrum_type": "ELNES"})
        self._check_run(d)

    def test_eels_wflow_abatom_by_symbol(self):
        wf_prim =get_wf_eels("O", self.structure, feff_input_set="ELNES",
                        edge="L1", user_tag_settings=self.user_tag_settings, use_primitive=True)
        wf = get_wf_eels("O", self.structure, feff_input_set="ELNES",
                        edge="L1", user_tag_settings=self.user_tag_settings, use_primitive=False)
        self.assertEqual(len(wf_prim.as_dict()["fws"]), 1)
        self.assertEqual(len(wf.as_dict()["fws"]), 2)

    def test_elnes_vs_exelfs(self):
        wf_elnes = get_wf_eels(self.absorbing_atom, self.structure, feff_input_set="ELNES",
                        edge="L1", user_tag_settings=self.user_tag_settings, use_primitive=True)
        wf_exelfs = get_wf_eels(self.absorbing_atom, self.structure, feff_input_set="EXELFS",
                        edge="L1", user_tag_settings=self.user_tag_settings, use_primitive=True)

        self.assertEqual(wf_elnes.as_dict()["fws"][0]["spec"]['_tasks'][0]['feff_input_set'][
                             '@class'], 'MPELNESSet')
        self.assertEqual(wf_exelfs.as_dict()["fws"][0]["spec"]['_tasks'][0]['feff_input_set'][
                             '@class'], 'MPEXELFSSet')

    def _check_run(self, d):
        run_dir = d["dir_name"]
        self.assertEqual(d["edge"], self.edge)
        self.assertEqual(d["absorbing_atom"], self.absorbing_atom)
        tags = Tags.from_file(os.path.join(run_dir, "feff.inp"))
        self.assertEqual(d["input_parameters"], tags.as_dict())


if __name__ == "__main__":
    unittest.main()
