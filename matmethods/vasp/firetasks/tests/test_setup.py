import os
from pymatgen import IStructure, Lattice
from pymatgen.io.vasp import Incar, Poscar, Potcar, Kpoints
from pymatgen.io.vasp.sets import MPVaspInputSet

import unittest

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

class TestSetup(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                                [1.9200989668, 3.3257101909, 0.00],
                                [0.00, -2.2171384943, 3.1355090603]])
        cls.struct_si = IStructure(lattice, ["Si"] * 2, coords)

        cls.ref_incar = Incar.from_file(os.path.join(module_dir, "reference_files", "setup_test", "INCAR"))
        cls.ref_poscar = Poscar.from_file(os.path.join(module_dir, "reference_files", "setup_test", "POSCAR"))
        cls.ref_potcar = Potcar.from_file(os.path.join(module_dir, "reference_files", "setup_test", "POTCAR"))
        cls.ref_kpoints = Kpoints.from_file(os.path.join(module_dir, "reference_files", "setup_test", "KPOINTS"))


    def tearDown(self):
        for x in ["INCAR", "POSCAR", "POTCAR", "KPOINTS"]:
            if os.path.exists(os.path.join(module_dir, x)):
                os.remove(os.path.join(module_dir, x))

    def _verify_files(self):
        self.assertEqual(Incar.from_file(os.path.join(module_dir, "INCAR")), self.ref_incar)
        self.assertEqual(str(Poscar.from_file(os.path.join(module_dir, "POSCAR"))), str(self.ref_poscar))
        self.assertEqual(Potcar.from_file(os.path.join(module_dir, "POTCAR")), self.ref_potcar)
        self.assertEqual(str(Kpoints.from_file(os.path.join(module_dir, "KPOINTS"))), str(self.ref_kpoints))

    def test_setup(self):
        try:
            mpvi = MPVaspInputSet()
            mpvi.write_input(self.struct_si, ".")
        except ValueError:
            import traceback
            traceback.print_exc()

            help_str = "This system is not set up to run VASP jobs. See further error tracebacks for help. Common fixes are (i) setting your VASP_PSP_DIR environment variable and (ii) making sure your VASP_PSP_DIR has the proper subdirs as outlined in PotcarSingle class of pymatgen, e.g. POT_GGA_PAW_PBE subdir."
            raise ValueError(help_str)

        self._verify_files()