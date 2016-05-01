import unittest

from matmethods.utils.loaders import get_wf_from_spec_dict
import os
from pymatgen import Structure
from pymatgen.util.testing import PymatgenTest

from monty.serialization import loadfn

class FuncTest(PymatgenTest):

    def setUp(self):
        self.structure = PymatgenTest.get_structure("Si")

    def test_get_wf_from_spec_dict(self):
        d = loadfn(os.path.join(os.path.abspath(os.path.dirname(__file__)), "spec.yaml"))
        wf = get_wf_from_spec_dict(self.structure, d)
        self.assertEqual(len(wf.fws), 4)
        for f in wf.fws:
            self.assertEqual(f.spec['_tasks'][-1]["db_file"], "db.json")

        self.assertEqual(sorted([len(v) for v in wf.links.values()]),
                         [0, 0, 1, 2])

if __name__ == '__main__':
    unittest.main()
