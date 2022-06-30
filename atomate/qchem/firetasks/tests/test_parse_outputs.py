import os
import shutil
import unittest

from pymatgen.core import Molecule
from pymatgen.io.qchem.inputs import QCInput

from atomate.qchem.firetasks.parse_outputs import QChemToDb
from atomate.utils.testing import AtomateTest

from monty.serialization import loadfn

__author__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"

module_dir = os.path.dirname(os.path.abspath(__file__))


class TestParseOutputQChem(AtomateTest):
    def test_parse_grad_good(self):
        my_calc_dir = os.path.join(module_dir, "..", "..", "test_files","parse_grad_good")
        ft = QChemToDb(calc_dir=my_calc_dir, parse_grad_file=True)
        ft.run_task({})
        task_doc = loadfn(os.path.join(my_calc_dir,"task.json"))
        self.assertEqual(task_doc["output"]["final_energy"], -274.6893362188)
        self.assertEqual(len(task_doc["output"]["precise_grad"]), 10)
        self.assertEqual(task_doc["output"]["precise_grad"][0], [0.0090906486788787, 0.016150932052898, 0.0054568671405536])
        self.assertEqual(task_doc["output"]["precise_grad"][-1], [0.0014495621906601, -0.0018570062958895, 0.0012478282193499])
        os.remove(os.path.join(my_calc_dir, "task.json"))

    def test_parse_grad_bad(self):
        my_calc_dir = os.path.join(module_dir, "..", "..", "test_files","parse_grad_bad")
        ft = QChemToDb(calc_dir=my_calc_dir, parse_grad_file=True)
        ft.run_task({})
        task_doc = loadfn(os.path.join(my_calc_dir,"task.json"))
        self.assertEqual(task_doc["warnings"]["grad_file_missing"], True)
        os.remove(os.path.join(my_calc_dir, "task.json"))

if __name__ == "__main__":
    unittest.main()
