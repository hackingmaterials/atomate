# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import shutil
import unittest

from pymatgen import Structure

from atomate.feff.workflows.xas import get_wf_exafs_paths

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "reference_files", "db_connections")


class TestEXAFSPaths(unittest.TestCase):

    def setUp(self):
        self.struct = Structure.from_file(os.path.join(module_dir, "..", "..", "test_files",
                                                       "feo_781777.json"))
        self.scratch_dir = os.path.join(module_dir, "scratch")
        if os.path.exists(self.scratch_dir):
            shutil.rmtree(self.scratch_dir)
        os.makedirs(self.scratch_dir)
        os.chdir(self.scratch_dir)
        wflow = get_wf_exafs_paths(0, self.struct, [[249, 0], [85, 0]], feff_cmd="feff",  db_file=None)
        self.wf_dict = wflow.as_dict()
        self.fw1_dict = self.wf_dict["fws"][0]
        self.fw2_dict = self.wf_dict["fws"][1]        
        if self.wf_dict["fws"][0]['name'] not in ['FeO-EXAFS-K-0']:
            self.fw1_dict, self.fw2_dict = self.fw2_dict, self.fw1_dict        

    def test_wflow_composition(self):
        self.assertEqual(len(self.wf_dict["fws"]), 2)
        ans = sorted(['FeO-EXAFS-K-0', 'FeO-EXAFS Paths'])
        self.assertEqual(ans, sorted([ft["name"] for ft in self.wf_dict["fws"]]))

    def test_feff_input_sets(self):
        ans_fis_fw1 = {'@class': 'MPEXAFSSet',
                       '@module': 'pymatgen.io.feff.sets',
                       'absorbing_atom': 0,
                       'edge': 'K',
                       'nkpts': 1000,
                       'radius': 10.0,
                       'user_tag_settings': {},
                       'structure': self.struct.as_dict()}
        ans_fis_fw2 = {'@class': 'MPEXAFSSet',
                       '@module': 'pymatgen.io.feff.sets',
                       'absorbing_atom': 0,
                       'edge': 'K',
                       'nkpts': 1000,
                       'radius': 10.0,
                       'user_tag_settings': {'CONTROL': '0 0 0 0 1 1', 'PRINT': '0 0 0 1 0 3'},
                       'structure': self.struct.as_dict()}
        fis_fw1 = self.fw1_dict["spec"]['_tasks'][0]['feff_input_set']
        fis_fw2 = self.fw2_dict["spec"]['_tasks'][1]['feff_input_set']
        self.assertDictEqual(fis_fw1, ans_fis_fw1)
        self.assertDictEqual(fis_fw2, ans_fis_fw2)

    def test_paths(self):
        paths = [[249, 0], [85, 0]]
        self.assertEqual(paths, self.fw2_dict["spec"]['_tasks'][2]['paths'])

    def tearDown(self):
        shutil.rmtree(self.scratch_dir)


if __name__ == "__main__":
    unittest.main()
