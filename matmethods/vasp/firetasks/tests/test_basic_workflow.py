import os

from monty.os.path import which

from fireworks.utilities.fw_serializers import load_object
from matmethods.vasp.examples.basic_vasp_workflows import get_basic_workflow
from matmethods.vasp.firetasks.write_inputs import WriteVaspFromIOSet, WriteVaspFromPMGObjects, ModifyIncar
from pymatgen import IStructure, Lattice
from pymatgen.io.vasp import Incar, Poscar, Potcar, Kpoints
from pymatgen.io.vasp.sets import MPVaspInputSet

import unittest

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
VASP_CMD = "vasp"

class BasicWorkflowTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if not which(VASP_CMD):
            raise unittest.SkipTest('VASP executable ("{}") cannot be found. SKIPPING {}.'.format(VASP_CMD, cls.__name__))

        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                                [1.9200989668, 3.3257101909, 0.00],
                                [0.00, -2.2171384943, 3.1355090603]])
        cls.struct_si = IStructure(lattice, ["Si"] * 2, coords)

        cls.ref_incar = Incar.from_file(os.path.join(module_dir, "reference_files", "INCAR"))
        cls.ref_poscar = Poscar.from_file(os.path.join(module_dir, "reference_files", "POSCAR"))
        cls.ref_potcar = Potcar.from_file(os.path.join(module_dir, "reference_files", "POTCAR"))
        cls.ref_kpoints = Kpoints.from_file(os.path.join(module_dir, "reference_files", "KPOINTS"))

    """
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

    """

    def test_relax_workflow(self):
        vis = MPVaspInputSet()
        structure = self.struct_si
        my_wf = get_basic_workflow(structure, vis, VASP_CMD)
        