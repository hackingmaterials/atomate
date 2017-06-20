# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
import copy
import ruamel.yaml as yaml

from fireworks.core.fworker import FWorker
from fireworks.core.rocket_launcher import rapidfire
from fireworks.utilities.fw_serializers import load_object
from atomate.vasp.firetasks.neb_tasks import TransferNEBTask, RunNEBVaspFake, \
    WriteNEBFromImages, WriteNEBFromEndpoints

from atomate.vasp.powerups import use_fake_vasp
from atomate.vasp.workflows.presets.core import wf_nudged_elastic_band
from atomate.utils.testing import AtomateTest

from pymatgen.core import Structure
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.vasp import Incar, Poscar, Potcar, Kpoints
from pymatgen.io.vasp.sets import MPRelaxSet

try:
    from pymatgen_diffusion.neb.io import get_endpoints_from_index

    pmgd = True
except ImportError:
    pmgd = False

__author__ = 'Hanmei Tang, Iek-Heng Chu'
__email__ = 'hat003@eng.ucsd.edu, ihchu@eng.ucsd.edu'

# /Users/hanmeiTang/repos/atomate/atomate/vasp/firetasks/tests  # todo: delete this line
module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class TestNEBTasks(PymatgenTest, AtomateTest):
    @classmethod
    def setUpClass(cls):
        # Structures used for test:
        cls.parent = PymatgenTest.get_structure("Li2O")

        cls.ref_incar = Incar.from_file(
            os.path.join(module_dir, "..", "..", "test_files", "neb_tasks_test", "INCAR"))
        cls.ref_poscar = Poscar.from_file(
            os.path.join(module_dir, "..", "..", "test_files", "neb_tasks_test", "POSCAR"))
        cls.ref_potcar = Potcar.from_file(
            os.path.join(module_dir, "..", "..", "test_files", "neb_tasks_test", "POTCAR"))
        cls.ref_kpoints = Kpoints.from_file(
            os.path.join(module_dir, "..", "..", "test_files", "neb_tasks_test", "KPOINTS"))
        cls.ref_incar_preserve = Incar.from_file(os.path.join(module_dir,
                                                              "..", "..", "test_files",
                                                              "neb_tasks_test", "INCAR"))

    def setUp(self):
        os.chdir(module_dir)

    def tearDown(self):
        for x in ["INCAR", "POSCAR", "POTCAR", "KPOINTS"]:
            if os.path.exists(os.path.join(module_dir, x)):
                os.remove(os.path.join(module_dir, x))

    #
    # def _verify_files(self, skip_kpoints=False, preserve_incar=False):
    #     if not preserve_incar:
    #         self.assertEqual(
    #             Incar.from_file(os.path.join(module_dir, "INCAR")),
    #             self.ref_incar)
    #         self.assertEqual(
    #             str(Poscar.from_file(os.path.join(module_dir, "POSCAR"))),
    #             str(self.ref_poscar))
    #         self.assertEqual((Potcar.from_file(os.path.join(module_dir,
    #                                                         "POTCAR"))).symbols,
    #                          self.ref_potcar.symbols)
    #         if not skip_kpoints:
    #             self.assertEqual(
    #                 str(Kpoints.from_file(
    #                     os.path.join(module_dir, "KPOINTS"))),
    #                 str(self.ref_kpoints))
    #     else:
    #         self.assertEqual(
    #             Incar.from_file(os.path.join(module_dir, "INCAR")),
    #             self.ref_incar_preserve)
    #
    # def test_ioset_explicit(self):
    #     ft = WriteVaspFromIOSet(dict(structure=self.struct_si, vasp_input_set=MPRelaxSet(self.struct_si, force_gamma=True)))
    #     ft = load_object(ft.to_dict())  # simulate database insertion
    #     ft.run_task({})
    #     self._verify_files()
    #
    # def test_ioset_implicit(self):
    #     ft = WriteVaspFromIOSet(
    #         dict(structure=self.struct_si, vasp_input_set="MPRelaxSet"))
    #     ft = load_object(ft.to_dict())  # simulate database insertion
    #     ft.run_task({})
    #     self._verify_files(skip_kpoints=True)
    #
    # def test_ioset_params(self):
    #     ft = WriteVaspFromIOSet(
    #         dict(structure=self.struct_si, vasp_input_set="MPRelaxSet",
    #              vasp_input_params={"user_incar_settings": {"ISMEAR": 1000}}))
    #     ft = load_object(ft.to_dict())  # simulate database insertion
    #     ft.run_task({})
    #     incar = Incar.from_file(os.path.join(module_dir, "INCAR"))
    #     self.assertEqual(incar["ISMEAR"], 1000)  # make sure override works
    #     incar['ISMEAR'] = -5  # switch back to default
    #     incar.write_file("INCAR")
    #     self._verify_files(skip_kpoints=True)
    #
    # def test_pmgobjects(self):
    #     mpvis = MPRelaxSet(self.struct_si, force_gamma=True)
    #     ft = WriteVaspFromPMGObjects({"incar": mpvis.incar,
    #                                   "poscar": mpvis.poscar,
    #                                   "kpoints": mpvis.kpoints,
    #                                   "potcar": mpvis.potcar})
    #     ft = load_object(ft.to_dict())  # simulate database insertion
    #     ft.run_task({})
    #     self._verify_files()
    #
    # def _get_processed_incar_dict(self, incar, poscar):
    #     incar_dict = incar.as_dict()
    #     comp = poscar.structure.composition
    #     elements = sorted([el for el in comp.elements if comp[el] > 0],
    #                       key=lambda e: e.X)
    #     most_electroneg = elements[-1].symbol
    #     for lda_param in ("LDAUL", "LDAUU", "LDAUJ"):
    #         if incar_dict.get(lda_param):
    #             vals = incar_dict[lda_param]
    #             if isinstance(vals, list):
    #                 incar_dict[lda_param] = {most_electroneg: {}}
    #                 for i, sym in enumerate(poscar.site_symbols):
    #                     incar_dict[lda_param][most_electroneg][sym] = vals[i]
    #     return incar_dict
    #
    # def test_modifyincar(self):
    #     # create an INCAR
    #     incar = self.ref_incar
    #     incar.write_file(os.path.join(module_dir, "INCAR"))
    #
    #     # modify and test
    #     ft = ModifyIncar(
    #         {"incar_update": {"ISMEAR": 1000}, "incar_multiply": {"ENCUT": 1.5},
    #          "incar_dictmod": {"_inc": {"ISPIN": -1}}})
    #     ft = load_object(ft.to_dict())  # simulate database insertion
    #     ft.run_task({})
    #
    #     incar_mod = Incar.from_file("INCAR")
    #     self.assertEqual(incar_mod['ISMEAR'], 1000)
    #     self.assertEqual(incar_mod['ENCUT'], 780)
    #     self.assertEqual(incar_mod['ISPIN'], 1)
    def test_transfer_neb_task(self):
        # TODO: Add unittest for class TransferNEBTask
        pass

    def test_run_neb_vasp_fake(self):
        # TODO: Add unittest for class RunNEBVaspFake
        pass

    def test_write_neb_from_images(self):
        # TODO: Add unittest for class WriteNEBFromImages
        pass

    def test_write_neb_from_endpoints(self):
        # TODO: Add unittest for class WriteNEBFromEndpoints
        pass


if __name__ == '__main__':
    unittest.main()
