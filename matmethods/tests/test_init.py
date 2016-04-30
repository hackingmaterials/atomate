import unittest

from matmethods import get_wf_from_spec_dict
import os
from pymatgen import Structure

from monty.serialization import loadfn

class FuncTest(unittest.TestCase):

    def setUp(self):
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = [[3.8401979337, 0.00, 0.00],
                   [1.9200989668, 3.3257101909, 0.00],
                   [0.00, -2.2171384943, 3.1355090603]]
        self.structure = Structure(lattice, ["Si"] * 2, coords)


    def test_get_wf_from_spec_dict(self):
        d = loadfn(os.path.join(os.path.abspath(os.path.dirname(__file__)), "spec.yaml"))
        wf = get_wf_from_spec_dict(self.structure, d)
        self.assertEqual(len(wf.fws), 4)
        for f in wf.fws:
            self.assertEqual(f.spec['_tasks'][-1]["db_file"], "db.json")
        self.assertEqual(wf.links,{-4: [], -1: [-2], -2: [-3, -4], -3: []})

if __name__ == '__main__':
    unittest.main()
