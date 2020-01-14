from pathlib import Path

from fireworks.utilities.fw_serializers import load_object

from atomate.vasp.firetasks.write_inputs import (
    WriteVaspFromIOSet,
    WriteVaspFromPMGObjects,
    ModifyPotcar,
    ModifyIncar,
)
from atomate.utils.testing import AtomateTest

from pymatgen.util.testing import PymatgenTest
from pymatgen.io.vasp import Incar, Poscar, Potcar, Kpoints
from pymatgen.io.vasp.sets import MPRelaxSet

__author__ = "Anubhav Jain, Kiran Mathew, Alex Ganose"
__email__ = "ajain@lbl.gov, kmathew@lbl.gov"

module_dir = Path(__file__).resolve().parent
test_files = module_dir / "../../test_files"
p_setup_test = test_files / "setup_test"
p_preserve_incar = test_files / "preserve_incar"


class TestWriteVasp(AtomateTest):
    @classmethod
    def setUpClass(cls):
        cls.struct_si = PymatgenTest.get_structure("Si")

        cls.ref_incar = Incar.from_file(p_setup_test / "INCAR")
        cls.ref_poscar = Poscar.from_file(p_setup_test / "POSCAR")
        cls.ref_potcar = Potcar.from_file(p_setup_test / "POTCAR")
        cls.ref_kpoints = Kpoints.from_file(p_setup_test / "KPOINTS")
        cls.ref_incar_preserve = Incar.from_file(p_preserve_incar / "INCAR")

    def setUp(self):
        super(TestWriteVasp, self).setUp(lpad=False)

    def tearDown(self):
        for x in ["INCAR", "POSCAR", "POTCAR", "KPOINTS", "POTCAR.spec"]:
            f = Path(self.scratch_dir) / x
            if f.exists():
                f.unlink()

    def _verify_files(
        self, skip_kpoints=False, preserve_incar=False, potcar_spec=False
    ):
        if not preserve_incar:
            self.assertEqual(Incar.from_file("INCAR"), self.ref_incar)

            poscar = Poscar.from_file("POSCAR")
            self.assertEqual(str(poscar), str(self.ref_poscar))

            if potcar_spec:
                symbols = Path("POTCAR.spec").read_text().split()
                self.assertEqual(symbols, self.ref_potcar.symbols)
            else:
                potcar = Potcar.from_file("POTCAR")
                self.assertEqual(potcar.symbols, self.ref_potcar.symbols)

            if not skip_kpoints:
                kpoints = Kpoints.from_file("KPOINTS")
                self.assertEqual(str(kpoints), str(self.ref_kpoints))
        else:
            self.assertEqual(Incar.from_file("INCAR"), self.ref_incar_preserve)

    def test_ioset_explicit(self):
        vis = MPRelaxSet(self.struct_si, force_gamma=True)
        ft = WriteVaspFromIOSet(structure=self.struct_si, vasp_input_set=vis)
        ft = load_object(ft.to_dict())  # simulate database insertion
        ft.run_task({})
        self._verify_files()

    def test_ioset_implicit(self):
        ft = WriteVaspFromIOSet(
            structure=self.struct_si, vasp_input_set="MPRelaxSet"
        )
        ft = load_object(ft.to_dict())  # simulate database insertion
        ft.run_task({})
        self._verify_files(skip_kpoints=True)

    def test_potcar_spec(self):
        ft = WriteVaspFromIOSet(
            structure=self.struct_si,
            vasp_input_set="MPRelaxSet",
            potcar_spec=True,
        )
        ft = load_object(ft.to_dict())  # simulate database insertion
        ft.run_task({})
        self._verify_files(potcar_spec=True)

    def test_ioset_params(self):
        ft = WriteVaspFromIOSet(
            structure=self.struct_si,
            vasp_input_set="MPRelaxSet",
            vasp_input_params={"user_incar_settings": {"ISMEAR": 1000}},
        )
        ft = load_object(ft.to_dict())  # simulate database insertion
        ft.run_task({})
        incar = Incar.from_file("INCAR")
        self.assertEqual(incar["ISMEAR"], 1000)  # make sure override works
        incar["ISMEAR"] = -5  # switch back to default
        incar.write_file("INCAR")
        self._verify_files(skip_kpoints=True)

    def test_pmgobjects(self):
        mpvis = MPRelaxSet(self.struct_si, force_gamma=True)
        ft = WriteVaspFromPMGObjects(
            incar=mpvis.incar,
            poscar=mpvis.poscar,
            kpoints=mpvis.kpoints,
            potcar=mpvis.potcar,
        )
        ft = load_object(ft.to_dict())  # simulate database insertion
        ft.run_task({})
        self._verify_files()

    def test_modify_incar(self):
        # create an INCAR
        incar = self.ref_incar
        incar.write_file("INCAR")

        # modify and test
        ft = ModifyIncar(
            incar_update={"ISMEAR": 1000},
            incar_multiply={"ENCUT": 1.5},
            incar_dictmod={"_inc": {"ISPIN": -1}},
        )
        ft = load_object(ft.to_dict())  # simulate database insertion
        ft.run_task({})

        incar_mod = Incar.from_file("INCAR")
        self.assertEqual(incar_mod["ISMEAR"], 1000)
        self.assertEqual(incar_mod["ENCUT"], 780)
        self.assertEqual(incar_mod["ISPIN"], 1)

    def test_modify_potcar(self):
        Potcar(["Si"]).write_file("POTCAR")
        potcar = Potcar.from_file("POTCAR")
        self.assertFalse("alt" in potcar[0].header)

        # modify/test
        ft = ModifyPotcar(potcar_symbols={"Si": "O"})
        ft = load_object(ft.to_dict())  # simulate database insertion
        ft.run_task({})

        new_potcar = Potcar.from_file("POTCAR")
        self.assertEqual(len(new_potcar), 1)
        self.assertEqual(new_potcar[0].symbol, "O")
