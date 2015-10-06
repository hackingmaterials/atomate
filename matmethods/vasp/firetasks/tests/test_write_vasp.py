import os
from fireworks.utilities.fw_serializers import load_object
from matmethods.vasp.firetasks.write_vasp import WriteVaspFromIOSet, WriteVaspFromPMGObjects, ModifyIncar
from pymatgen import IStructure, Lattice
from pymatgen.io.vasp import Incar, Poscar, Potcar, Kpoints
from pymatgen.io.vasp.sets import MPVaspInputSet

import unittest

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

class WriteVaspTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                                [1.9200989668, 3.3257101909, 0.00],
                                [0.00, -2.2171384943, 3.1355090603]])
        cls.struct = IStructure(lattice, ["Si"] * 2, coords)

        cls.ref_incar = Incar.from_file(os.path.join(module_dir, "reference_files", "INCAR"))
        cls.ref_poscar = Poscar.from_file(os.path.join(module_dir, "reference_files", "POSCAR"))
        cls.ref_potcar = Potcar.from_file(os.path.join(module_dir, "reference_files", "POTCAR"))
        cls.ref_kpoints = Kpoints.from_file(os.path.join(module_dir, "reference_files", "KPOINTS"))

    @classmethod
    def tearDownClass(cls):
        pass

    def setUp(self):
        pass

    def tearDown(self):
        for x in ["INCAR", "POSCAR", "POTCAR", "KPOINTS"]:
            if os.path.exists(os.path.join(module_dir, x)):
                os.remove(os.path.join(module_dir, x))

    def _verify_files(self):
        self.assertEqual(Incar.from_file(os.path.join(module_dir, "INCAR")), self.ref_incar)
        self.assertEqual(str(Poscar.from_file(os.path.join(module_dir, "POSCAR"))), str(self.ref_poscar))
        self.assertEqual(Potcar.from_file(os.path.join(module_dir, "POTCAR")), self.ref_potcar)
        self.assertEqual(str(Kpoints.from_file(os.path.join(module_dir, "KPOINTS"))), str(self.ref_kpoints))

    def test_ioset_explicit(self):
        ft = WriteVaspFromIOSet(dict(structure=self.struct, vasp_input_set=MPVaspInputSet()))
        ft = load_object(ft.to_dict())  # simulate database insertion
        ft.run_task({})
        self._verify_files()

    def test_ioset_implicit(self):
        ft = WriteVaspFromIOSet(dict(structure=self.struct, vasp_input_set="MPVaspInputSet"))
        ft = load_object(ft.to_dict())  # simulate database insertion
        ft.run_task({})
        self._verify_files()

    def test_ioset_params(self):
        ft = WriteVaspFromIOSet(dict(structure=self.struct, vasp_input_set="MPVaspInputSet",
                                     vasp_input_params={"user_incar_settings": {"ISMEAR": 1000}}))
        ft = load_object(ft.to_dict())  # simulate database insertion
        ft.run_task({})
        incar = Incar.from_file(os.path.join(module_dir, "INCAR"))
        self.assertEqual(incar["ISMEAR"], 1000)  # make sure override works
        incar['ISMEAR'] = -5  # switch back to default
        incar.write_file("INCAR")
        self._verify_files()

    def test_pmgobjects(self):
        mpvis = MPVaspInputSet()
        ft = WriteVaspFromPMGObjects({"incar": mpvis.get_incar(self.struct),
                                      "poscar": mpvis.get_poscar(self.struct),
                                      "kpoints": mpvis.get_kpoints(self.struct),
                                      "potcar": mpvis.get_potcar(self.struct)})
        ft = load_object(ft.to_dict())  # simulate database insertion
        ft.run_task({})
        self._verify_files()

    def test_modifyincar(self):
        # create an INCAR
        incar = self.ref_incar
        incar.write_file(os.path.join(module_dir, "INCAR"))

        # modify and test
        ft = ModifyIncar({"key_update": {"ISMEAR": 1000}, "key_multiply": {"ENCUT": 1.5}, "key_dictmod": {"_inc": {"ISPIN": -1}}})
        ft = load_object(ft.to_dict())  # simulate database insertion
        ft.run_task({})

        incar_mod = Incar.from_file("INCAR")
        self.assertEqual(incar_mod['ISMEAR'], 1000)
        self.assertEqual(incar_mod['ENCUT'], 780)
        self.assertEqual(incar_mod['ISPIN'], 1)



if __name__ == '__main__':
    unittest.main()
