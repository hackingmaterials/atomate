# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

from pymatgen import IStructure, Lattice
from pymatgen.io.vasp import Incar, Poscar, Potcar, Kpoints
from pymatgen.io.vasp.sets import MPRelaxSet

__author__ = 'Anubhav Jain, Kiran Mathew'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class TestSetup(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        cls.struct_si = IStructure(lattice, ["Si"] * 2, coords)

        cls.ref_incar = Incar.from_file(
            os.path.join(module_dir, "..", "test_files", "setup_test", "INCAR"))
        cls.ref_poscar = Poscar.from_file(
            os.path.join(module_dir, "..", "test_files", "setup_test", "POSCAR"))
        cls.ref_potcar = Potcar.from_file(
            os.path.join(module_dir, "..", "test_files", "setup_test", "POTCAR"))
        cls.ref_kpoints = Kpoints.from_file(
            os.path.join(module_dir, "..", "test_files", "setup_test", "KPOINTS"))

    def setUp(self):
        os.chdir(module_dir)

    def tearDown(self):
        for x in ["INCAR", "POSCAR", "POTCAR", "KPOINTS"]:
            if os.path.exists(os.path.join(module_dir, x)):
                os.remove(os.path.join(module_dir, x))

    def _verify_files(self):
        self.assertEqual(Incar.from_file(os.path.join(module_dir, "INCAR")), self.ref_incar)
        self.assertEqual(str(Poscar.from_file(os.path.join(module_dir, "POSCAR"))),
                         str(self.ref_poscar))
        self.assertEqual(Potcar.from_file(os.path.join(module_dir, "POTCAR")).symbols,
                         self.ref_potcar.symbols)
        self.assertEqual(str(Kpoints.from_file(os.path.join(module_dir, "KPOINTS"))),
                         str(self.ref_kpoints))

    def test_setup(self):
        try:
            vi = MPRelaxSet(self.struct_si, force_gamma=True)
            vi.write_input(".")
        except ValueError:
            import traceback
            traceback.print_exc()

            help_str = "This system is not set up to run VASP jobs. See further error tracebacks " \
                       "for help. Try making sure your PMG_VASP_PSP_DIR has the proper subdirs as " \
                       "outlined in PotcarSingle class of pymatgen, e.g. POT_GGA_PAW_PBE subdir."
            raise ValueError(help_str)

        self._verify_files()


if __name__ == "__main__":
    unittest.main()
