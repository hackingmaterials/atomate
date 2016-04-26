# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

from fireworks.utilities.fw_serializers import load_object

from matmethods.vasp.firetasks.write_inputs import WriteVaspFromIOSet, WriteVaspFromPMGObjects, ModifyIncar
from matmethods.vasp.input_sets import StructureOptimizationVaspInputSet, StaticVaspInputSet, write_with_preserved_incar

from pymatgen import IStructure, Lattice, Structure
from pymatgen.io.vasp import Incar, Poscar, Potcar, Kpoints

__author__ = 'Anubhav Jain, Kiran Mathew'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class TestWriteVasp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not os.environ.get("VASP_PSP_DIR"):
            raise unittest.SkipTest(
                'This system is not set up to run VASP jobs. '
                'Please set your VASP_PSP_DIR environment variable.')

        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        cls.struct_si = IStructure(lattice, ["Si"] * 2, coords)

        cls.ref_incar = Incar.from_file(
            os.path.join(module_dir, "reference_files", "setup_test", "INCAR"))
        cls.ref_poscar = Poscar.from_file(
            os.path.join(module_dir, "reference_files", "setup_test",
                         "POSCAR"))
        cls.ref_potcar = Potcar.from_file(
            os.path.join(module_dir, "reference_files", "setup_test",
                         "POTCAR"))
        cls.ref_kpoints = Kpoints.from_file(
            os.path.join(module_dir, "reference_files", "setup_test",
                         "KPOINTS"))
        cls.ref_incar_preserve = Incar.from_file(os.path.join(module_dir,
                                                              "reference_files",
                                                              "preserve_incar",
                                                              "INCAR"))
        cls.ref_incar_preserve.update(StaticVaspInputSet.DEFAULT_SETTINGS)

    def setUp(self):
        os.chdir(module_dir)

    def tearDown(self):
        for x in ["INCAR", "POSCAR", "POTCAR", "KPOINTS"]:
            if os.path.exists(os.path.join(module_dir, x)):
                os.remove(os.path.join(module_dir, x))

    def _verify_files(self, skip_kpoints=False, preserve_incar=False):
        if not preserve_incar:
            self.assertEqual(
                Incar.from_file(os.path.join(module_dir, "INCAR")),
                self.ref_incar)
            self.assertEqual(
                str(Poscar.from_file(os.path.join(module_dir, "POSCAR"))),
                str(self.ref_poscar))
            self.assertEqual((Potcar.from_file(os.path.join(module_dir,
                                                            "POTCAR"))).symbols,
                             self.ref_potcar.symbols)
            if not skip_kpoints:
                self.assertEqual(
                    str(Kpoints.from_file(
                        os.path.join(module_dir, "KPOINTS"))),
                    str(self.ref_kpoints))
        else:
            self.assertEqual(
                Incar.from_file(os.path.join(module_dir, "INCAR")),
                self.ref_incar_preserve)

    def test_ioset_explicit(self):
        ft = WriteVaspFromIOSet(dict(structure=self.struct_si,
                                     vasp_input_set=StructureOptimizationVaspInputSet()))
        ft = load_object(ft.to_dict())  # simulate database insertion
        ft.run_task({})
        self._verify_files()

    def test_ioset_implicit(self):
        ft = WriteVaspFromIOSet(
            dict(structure=self.struct_si, vasp_input_set="MPVaspInputSet"))
        ft = load_object(ft.to_dict())  # simulate database insertion
        ft.run_task({})
        self._verify_files(skip_kpoints=True)

    def test_ioset_params(self):
        ft = WriteVaspFromIOSet(
            dict(structure=self.struct_si, vasp_input_set="MPVaspInputSet",
                 vasp_input_params={"user_incar_settings": {"ISMEAR": 1000}}))
        ft = load_object(ft.to_dict())  # simulate database insertion
        ft.run_task({})
        incar = Incar.from_file(os.path.join(module_dir, "INCAR"))
        self.assertEqual(incar["ISMEAR"], 1000)  # make sure override works
        incar['ISMEAR'] = -5  # switch back to default
        incar.write_file("INCAR")
        self._verify_files(skip_kpoints=True)

    def test_pmgobjects(self):
        mpvis = StructureOptimizationVaspInputSet()
        ft = WriteVaspFromPMGObjects({"incar": mpvis.get_incar(self.struct_si),
                                      "poscar": mpvis.get_poscar(
                                          self.struct_si),
                                      "kpoints": mpvis.get_kpoints(
                                          self.struct_si),
                                      "potcar": mpvis.get_potcar(
                                          self.struct_si)})
        ft = load_object(ft.to_dict())  # simulate database insertion
        ft.run_task({})
        self._verify_files()

    def test_preserve_incar(self):
        prev_structure = Structure.from_file(os.path.join(module_dir,
                                                          "reference_files",
                                                          "preserve_incar",
                                                          "POSCAR_inverted"))
        prev_structure_decorated = prev_structure.copy()
        prev_structure_decorated.add_site_property("magmom",
                                                   [0.0, 0.0, 2.0, -2.0])
        # new structure is sorted and is a supercell of the previous one
        new_structure = prev_structure_decorated.copy()
        new_structure.sort()
        new_structure.make_supercell([2, 1, 1])
        prev_incar = Incar.from_file(os.path.join(module_dir,
                                                  "reference_files",
                                                  "preserve_incar",
                                                  "INCAR_inverted"))
        # overide the ldau params read from the default yaml inputset
        # get_incar method expects ldau params in {"most_electroneg":{
        # "symbol": value}} format and incar.as_dict() yields the ldau
        # params as list
        config_dict = {
            "INCAR": self._get_processed_incar_dict(prev_incar,
                                                    Poscar(prev_structure))}
        vis = StaticVaspInputSet()
        vis.incar_settings = config_dict['INCAR']
        vis.incar_settings.update(vis.DEFAULT_SETTINGS)
        write_with_preserved_incar(vis, new_structure, os.path.join(module_dir,
                                                               "reference_files",
                                                               "preserve_incar"))
        self._verify_files(preserve_incar=True)
        # test MAGMOM with LSORBIT
        config_dict_override = {"INCAR": {"LSORBIT": "T"}}
        vis = StaticVaspInputSet()
        vis.incar_settings = config_dict['INCAR']
        vis.incar_settings.update(vis.DEFAULT_SETTINGS)
        vis.incar_settings.update({})
        # the incar in prev_dir has MAGMOM   = 2*2.0 2*-2.0 4*0
        write_with_preserved_incar(vis, new_structure,
                                   os.path.join(module_dir,"reference_files","preserve_incar"),
                                   config_dict_override=config_dict_override, output_dir='.')
        incar = Incar.from_file(os.path.join(module_dir, 'INCAR'))
        self.assertTrue(incar["LSORBIT"])
        self.assertEqual(incar["MAGMOM"], [[0, 0, 2.0],
                                           [0, 0, 2.0],
                                           [0, 0, -2.0],
                                           [0, 0, -2.0],
                                           [0, 0, 0.0],
                                           [0, 0, 0.0],
                                           [0, 0, 0.0],
                                           [0, 0, 0.0]])

    def _get_processed_incar_dict(self, incar, poscar):
        incar_dict = incar.as_dict()
        comp = poscar.structure.composition
        elements = sorted([el for el in comp.elements if comp[el] > 0],
                          key=lambda e: e.X)
        most_electroneg = elements[-1].symbol
        for lda_param in ("LDAUL", "LDAUU", "LDAUJ"):
            if incar_dict.get(lda_param):
                vals = incar_dict[lda_param]
                if isinstance(vals, list):
                    incar_dict[lda_param] = {most_electroneg: {}}
                    for i, sym in enumerate(poscar.site_symbols):
                        incar_dict[lda_param][most_electroneg][sym] = vals[i]
        return incar_dict

    def test_modifyincar(self):
        # create an INCAR
        incar = self.ref_incar
        incar.write_file(os.path.join(module_dir, "INCAR"))

        # modify and test
        ft = ModifyIncar(
            {"incar_update": {"ISMEAR": 1000}, "incar_multiply": {"ENCUT": 1.5},
             "incar_dictmod": {"_inc": {"ISPIN": -1}}})
        ft = load_object(ft.to_dict())  # simulate database insertion
        ft.run_task({})

        incar_mod = Incar.from_file("INCAR")
        self.assertEqual(incar_mod['ISMEAR'], 1000)
        self.assertEqual(incar_mod['ENCUT'], 780)
        self.assertEqual(incar_mod['ISPIN'], 1)


if __name__ == '__main__':
    unittest.main()
