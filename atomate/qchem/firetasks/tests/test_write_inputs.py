# coding: utf-8


import os
import unittest
import shutil

from atomate.qchem.firetasks.write_inputs import WriteInputFromIOSet, WriteInput, WriteCustomInput
from atomate.utils.testing import AtomateTest
from pymatgen.core import Molecule
from pymatgen.io.qchem.inputs import QCInput

__author__ = "Brandon Wood"
__email__ = "b.wood@berkeley.edu"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class TestWriteInputQChem(AtomateTest):
    @classmethod
    def setUpClass(cls):

        co_species = ["C", "O"]
        co_coords = [[0.0, 0.0, 0.0], [1.3, 0.0, 0.0]]
        cls.co_mol = Molecule(co_species, co_coords)
        cls.co_opt_ref_in = QCInput.from_file(
            os.path.join(module_dir, "..", "..", "test_files", "co_qc.in"))
        cls.opt_mol_ref_in = QCInput.from_file(
            os.path.join(module_dir, "..", "..", "test_files", "to_opt.qin"))
        cls.opt_mol = cls.opt_mol_ref_in.molecule
        cls.opt_mol_pcm_ref_in = QCInput.from_file(
            os.path.join(module_dir, "..", "..", "test_files",
                         "to_opt_pcm.qin"))
        cls.opt_mol_smd_ref_in = QCInput.from_file(
            os.path.join(module_dir, "..", "..", "test_files",
                         "to_opt_smd.qin"))

    def setUp(self, lpad=False):
        super(TestWriteInputQChem, self).setUp(lpad=False)

    def tearDown(self):
        shutil.rmtree(self.scratch_dir)
        for x in ["mol.qin"]:
            if os.path.exists(os.path.join(module_dir, x)):
                os.remove(os.path.join(module_dir, x))

    def test_write_input_from_io_set(self):
        ft = WriteInputFromIOSet(
            molecule=self.co_mol, qchem_input_set="OptSet")
        ft.run_task({})
        test_dict = QCInput.from_file("mol.qin").as_dict()
        for k, v in self.co_opt_ref_in.as_dict().items():
            self.assertEqual(v, test_dict[k])

    def test_write_input_from_io_set_diff_mol(self):
        ft = WriteInputFromIOSet(
            molecule=self.opt_mol, qchem_input_set="OptSet")
        ft.run_task({})
        test_dict = QCInput.from_file("mol.qin").as_dict()
        for k, v in self.opt_mol_ref_in.as_dict().items():
            self.assertEqual(v, test_dict[k])

    def test_write_input_from_io_set_diff_mol_pcm(self):
        ft = WriteInputFromIOSet(
            molecule=self.opt_mol,
            qchem_input_set="OptSet",
            qchem_input_params={"pcm_dielectric": 10.0})
        ft.run_task({})
        test_dict = QCInput.from_file("mol.qin").as_dict()
        for k, v in self.opt_mol_pcm_ref_in.as_dict().items():
            self.assertEqual(v, test_dict[k])

    def test_write_input_from_io_custom_smd(self):
        ft = WriteInputFromIOSet(
            molecule=self.opt_mol,
            qchem_input_set="OptSet",
            qchem_input_params={"smd_solvent": "custom", "custom_smd": "90.00,1.415,0.00,0.735,20.2,0.00,0.00"})
        ft.run_task({})
        test_dict = QCInput.from_file("mol.qin").as_dict()
        for k, v in self.opt_mol_smd_ref_in.as_dict().items():
            self.assertEqual(v, test_dict[k])
        with open("solvent_data") as sd:
            lines = sd.readlines()
            self.assertEqual(lines[0],"90.00,1.415,0.00,0.735,20.2,0.00,0.00")
        os.remove("solvent_data")

    def test_write_input_from_io_set_write_dir(self):
        ft = WriteInputFromIOSet(
            molecule=self.co_mol,
            qchem_input_set="OptSet",
            write_to_dir=module_dir)
        ft.run_task({})
        test_dict = QCInput.from_file(os.path.join(module_dir,
                                                   "mol.qin")).as_dict()
        for k, v in self.co_opt_ref_in.as_dict().items():
            self.assertEqual(v, test_dict[k])

    def test_write_input(self):
        mol = self.co_mol
        rem = {
            "job_type": "opt",
            "basis": "def2-tzvppd",
            "max_scf_cycles": 200,
            "method": "wB97xd",
            "geom_opt_max_cycles": 200,
            "gen_scfman": True,
            "scf_algorithm": "diis",
            "xc_grid": 3,
            "sym_ignore": True,
            "symmetry": False,
            "resp_charges": True
        }
        qc_input = QCInput(mol, rem)
        ft = WriteInput(qc_input=qc_input)
        ft.run_task({})
        test_dict = QCInput.from_file("mol.qin").as_dict()
        for k, v in self.co_opt_ref_in.as_dict().items():
            self.assertEqual(v, test_dict[k])

    def test_write_custom_input(self):
        mol = self.co_mol
        rem = {
            "job_type": "opt",
            "basis": "def2-tzvppd",
            "max_scf_cycles": 200,
            "method": "wB97xd",
            "geom_opt_max_cycles": 200,
            "gen_scfman": True,
            "scf_algorithm": "diis",
            "xc_grid": 3,
            "sym_ignore": True,
            "symmetry": False,
            "resp_charges": True
        }
        ft = WriteCustomInput(molecule=mol, rem=rem)
        ft.run_task({})
        test_dict = QCInput.from_file("mol.qin").as_dict()
        for k, v in self.co_opt_ref_in.as_dict().items():
            self.assertEqual(v, test_dict[k])


if __name__ == "__main__":
    unittest.main()
