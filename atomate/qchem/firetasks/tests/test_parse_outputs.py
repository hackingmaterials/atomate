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


class TestParseOutputQChem_grads(AtomateTest):
    def setUp(self, lpad=False):
        super().setUp(lpad=False)

    def tearDown(self):
        pass

    def test_parse_grad_good(self):
        my_calc_dir = os.path.join(module_dir, "..", "..", "test_files","parse_grad_good")
        ft = QChemToDb(calc_dir=my_calc_dir, parse_grad_file=True)
        ft.run_task({})
        task_doc = loadfn(os.path.join(my_calc_dir,"task.json"))
        self.assertEqual(task_doc["output"]["final_energy"], -274.6893362188)
        self.assertEqual(len(task_doc["output"]["precise_gradients"]), 10)
        self.assertEqual(task_doc["output"]["precise_gradients"][0], [0.0090906486788787, 0.016150932052898, 0.0054568671405536])
        self.assertEqual(task_doc["output"]["precise_gradients"][-1], [0.0014495621906601, -0.0018570062958895, 0.0012478282193499])
        os.remove(os.path.join(my_calc_dir, "task.json"))

    def test_parse_grad_131(self):
        my_calc_dir = os.path.join(module_dir, "..", "..", "test_files","tmqm_grad_pcm")
        ft = QChemToDb(calc_dir=my_calc_dir, parse_grad_file=True)
        ft.run_task({})
        task_doc = loadfn(os.path.join(my_calc_dir,"task.json"))
        self.assertEqual(task_doc["output"]["final_energy"], -2791.8404057999)
        self.assertEqual(len(task_doc["output"]["precise_gradients"]), 25)
        self.assertEqual(task_doc["output"]["precise_gradients"][0], [-2.7425178677332305e-05, 1.8017443772144412e-06, -2.3689773769176685e-06])
        self.assertEqual(task_doc["output"]["precise_gradients"][-1], [0.0028753270363098644, -0.000392640066359285, 0.004405091870637312])
        os.remove(os.path.join(my_calc_dir, "task.json"))

    def test_parse_grad_bad(self):
        my_calc_dir = os.path.join(module_dir, "..", "..", "test_files","parse_grad_bad")
        ft = QChemToDb(calc_dir=my_calc_dir, parse_grad_file=True)
        ft.run_task({})
        task_doc = loadfn(os.path.join(my_calc_dir,"task.json"))
        self.assertEqual(task_doc["warnings"]["grad_file_missing"], True)
        os.remove(os.path.join(my_calc_dir, "task.json"))


class TestParseOutputQChem_hess(AtomateTest):
    def setUp(self, lpad=False):
        os.makedirs(os.path.join(module_dir, "..", "..", "test_files", "freq_save_hess", "scratch"))
        shutil.copyfile(
            os.path.join(module_dir, "..", "..", "test_files", "freq_save_hess", "BUP_scratch", "132.0"),
            os.path.join(module_dir, "..", "..", "test_files", "freq_save_hess", "scratch", "132.0"),
        )
        shutil.copyfile(
            os.path.join(module_dir, "..", "..", "test_files", "freq_save_hess", "BUP_scratch", "HESS"),
            os.path.join(module_dir, "..", "..", "test_files", "freq_save_hess", "scratch", "HESS"),
        )
        super().setUp(lpad=False)

    def test_parse_hess(self):
        my_calc_dir = os.path.join(module_dir, "..", "..", "test_files","freq_save_hess")
        ft = QChemToDb(calc_dir=my_calc_dir, parse_hess_file=True)
        ft.run_task({})
        task_doc = loadfn(os.path.join(my_calc_dir,"task.json"))
        self.assertEqual(task_doc["output"]["final_energy"], -151.3244603665)
        self.assertEqual(task_doc["output"]["hess_data"]["scratch/132.0"][0],0.12636293260949633)
        self.assertEqual(task_doc["output"]["hess_data"]["scratch/132.0"][-2],-0.2025032138024329)
        self.assertEqual(task_doc["output"]["hess_data"]["scratch/HESS"][-2], ' -0.175476533300377      -0.202503213802433       0.205623571433770     \n')
        os.remove(os.path.join(my_calc_dir, "task.json"))

    def tearDown(self):
        pass

if __name__ == "__main__":
    unittest.main()
