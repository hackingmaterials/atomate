__author__ = "Janine George, Guido Petretto"
__email__ = "janine.george@uclouvain.be, guido.petretto@uclouvain.be"

import json
import os
import unittest

from atomate.utils.testing import AtomateTest
from atomate.utils.testing import DB_DIR
from atomate.vasp.powerups import use_fake_lobster
from atomate.vasp.powerups import use_fake_vasp, use_custodian
from atomate.vasp.workflows.base.lobster import get_wf_lobster, \
    get_wf_lobster_test_basis
from fireworks.core.rocket_launcher import rapidfire
from pymatgen.core.structure import Structure
from pymatgen.util.testing import PymatgenTest
from atomate.vasp.powerups import use_potcar_spec

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

VASP_CMD = None
#LOBSTER_CMD = None
DEBUG = False
DB_FILE = os.path.join(DB_DIR, "db.json")
refs_dirs_si_vasp = {"static": os.path.join(module_dir, "../../test_files/lobster/si_vasp_lobster/vasp/")}
refs_dirs_si_lobster = {
    "lobster_calculation": os.path.join(module_dir, "../../test_files/lobster/si_vasp_lobster/lobster")}
refs_dirs_complex_vasp = {
    "static": os.path.join(module_dir, "../../test_files/lobster/complex_vasp_lobster/vasp")
}

refs_dirs_complex_lobster = {
    "lobster_calculation_0": os.path.join(module_dir,
                                          "../../test_files/lobster/complex_vasp_lobster/lobster_0"),
    "lobster_calculation_1": os.path.join(module_dir, "../../test_files/lobster/complex_vasp_lobster/lobster_1")
}


class TestWFLobster(AtomateTest):
    def setUp(self):
        super(TestWFLobster, self).setUp(lpad=True)
        self.struct_si = Structure.from_file(os.path.join(refs_dirs_si_vasp["static"], "inputs", "POSCAR"))

    def _check_run(self, d, mode, database=True):
        if mode not in ["lobsternormal"]:
            raise ValueError("Invalid mode!")

        self.assertEqual(d["formula_pretty"], "Si")
        self.assertEqual(d["formula_anonymous"], "A")
        self.assertEqual(d["nelements"], 1)
        self.assertEqual(d["state"], "successful")

        if mode in ["lobsternormal", "lobsternormal_delete_wavecars"]:
            self.assertListEqual(d["output"]["chargespilling"], [0.0147, 0.0147])
            self.assertListEqual(d["output"]["elements"], ['Si'])
            self.assertListEqual(d["output"]["basistype"], ['pbeVaspFit2015'])
            self.assertListEqual(d["output"]["basisfunctions"], [['3s', '3p_y', '3p_z', '3p_x']])
            self.assertTrue(d["output"]["hasDOSCAR"])
            if database:
                self.assertNotIn("cohpcar_id", d)
                self.assertIn("icooplist_id", d)
            launch_dir = d["dir_name"]
            wavecar_present = os.path.isfile(os.path.join(launch_dir, "WAVECAR.gz")) or \
                              os.path.isfile(os.path.join(launch_dir, "WAVECAR"))
            if mode == "lobsternormal_delete_wavecars":
                self.assertFalse(wavecar_present)
            else:
                self.assertTrue(wavecar_present)

    def _single_vasp_lobster(self, delete_wavecars=False,
                             user_supplied_basis=None, fake=True):
        # add the workflow
        structure = self.struct_si
        my_wf = get_wf_lobster(structure=structure, c={"vasp_cmd": VASP_CMD, "DB_FILE": None},
                               user_kpoints_settings={"grid_density": 100}, delete_all_wavecars=delete_wavecars,
                               user_supplied_basis=user_supplied_basis)
        if fake:
            my_wf = use_fake_vasp(my_wf, refs_dirs_si_vasp)
            my_wf = use_fake_lobster(my_wf, refs_dirs_si_lobster)
        else:
            my_wf = use_custodian(my_wf)

        self.lp.add_wf(my_wf)

        # run the workflow
        rapidfire(self.lp)

        fw = self.lp.get_fw_by_id(fw_id=1)
        with open(os.path.join(fw.launches[-1].launch_dir, "task_lobster.json")) as f:
            d = json.load(f)
        self._check_run(d, mode="lobsternormal", database=False)

        wf = self.lp.get_wf_by_fw_id(1)
        self.assertTrue(all([s == 'COMPLETED' for s in wf.fw_states.values()]))

    def _single_lobster_db_insertion(self, delete_wavecars=False,
                                     user_supplied_basis=None, fake=True):

        structure = self.struct_si
        my_wf = get_wf_lobster(structure=structure, c={"vasp_cmd": VASP_CMD, "DB_FILE": DB_FILE},
                               user_kpoints_settings={"grid_density": 100}, delete_all_wavecars=delete_wavecars,
                               user_supplied_basis=user_supplied_basis,
                               additional_outputs=["ICOOPLIST.lobster"])
        if fake:
            my_wf = use_fake_vasp(my_wf, refs_dirs_si_vasp)
            my_wf = use_fake_lobster(my_wf, refs_dirs_si_lobster)
        else:
            my_wf = use_custodian(my_wf)

        self.lp.add_wf(my_wf)

        # run the workflow
        rapidfire(self.lp)
        d = self.get_task_collection(coll_name="lobster").find_one()
        if delete_wavecars:
            self._check_run(d, mode="lobsternormal_delete_wavecars")
        else:
            self._check_run(d, mode="lobsternormal")

        wf = self.lp.get_wf_by_fw_id(1)
        self.assertTrue(all([s == 'COMPLETED' for s in wf.fw_states.values()]))

    def test_single_vasp_lobster(self):
        self._single_vasp_lobster(fake=True)

    def test_single_lobster_db_insertion(self):
        self._single_lobster_db_insertion(fake=True)

    # def test_single_vasp_lobster_integration(self):
    #     # integration test
    #     if VASP_CMD and LOBSTER_CMD:
    #        self._single_vasp_lobster(fake=False)



class TestWFLobsterTestBasis(AtomateTest):
    """
    Test for get_wf_lobster_test_basis
    """

    def setUp(self):
        super(TestWFLobsterTestBasis, self).setUp(lpad=True)
        self.struct_mp = Structure.from_file(os.path.join(refs_dirs_complex_vasp["static"], "inputs", "POSCAR"))

    def _check_run(self, d, mode):
        if mode not in ["lobsternormal"]:
            raise ValueError("Invalid mode!")

        self.assertEqual(d["formula_pretty"], "CdF2")
        self.assertEqual(d["formula_anonymous"], "AB2")
        self.assertEqual(d["nelements"], 2)
        self.assertEqual(d["state"], "successful")

        if mode in ["lobsternormal"]:
            self.assertListEqual(d["output"]["chargespilling"], [0.0027, 0.0027])
            self.assertListEqual(d["output"]["elements"], ['F', 'Cd'])
            self.assertListEqual(d["output"]["basistype"], ['pbeVaspFit2015', 'pbeVaspFit2015'])
            self.assertListEqual(d["output"]["basisfunctions"], [['2s', '2p_y', '2p_z', '2p_x'],
                                                                 ['5s', '5p_y', '5p_z', '5p_x', '4d_xy', '4d_yz',
                                                                  '4d_z^2', '4d_xz', '4d_x^2-y^2']])
            self.assertTrue(d["output"]["hasDOSCAR"])

    def _single_vasp_lobster(self, fake=True):
        # add the workflow

        structure = self.struct_mp
        my_wf = get_wf_lobster_test_basis(structure=structure, c={"vasp_cmd": VASP_CMD, "DB_FILE": None},
                                          user_kpoints_settings={"grid_density": 100},
                                          delete_all_wavecars=False)
        if fake:
            my_wf = use_fake_vasp(my_wf, refs_dirs_complex_vasp)
            my_wf = use_fake_lobster(my_wf, ref_dirs=refs_dirs_complex_lobster)
        else:
            my_wf = use_custodian(my_wf)

        self.lp.add_wf(my_wf)

        # run the workflow
        rapidfire(self.lp)

        fw = self.lp.get_fw_by_id(fw_id=1)

        with open(os.path.join(fw.launches[-1].launch_dir, "task_lobster.json")) as f:
            d = json.load(f)
        self._check_run(d, mode="lobsternormal")

        wf = self.lp.get_wf_by_fw_id(1)
        self.assertTrue(all([s == 'COMPLETED' for s in wf.fw_states.values()]))

    def _single_lobster_db_insertion(self, fake=True):
        structure = self.struct_mp
        my_wf = get_wf_lobster_test_basis(structure=structure, c={"vasp_cmd": VASP_CMD, "DB_FILE": DB_FILE},
                                          user_kpoints_settings={"grid_density": 100},
                                          delete_all_wavecars=False)

        if fake:
            my_wf = use_fake_vasp(my_wf, ref_dirs=refs_dirs_complex_vasp)
            my_wf = use_fake_lobster(my_wf, ref_dirs=refs_dirs_complex_lobster)
        else:
            my_wf = use_custodian(my_wf)

        self.lp.add_wf(my_wf)

        # run the workflow
        rapidfire(self.lp)
        d = self.get_task_collection(coll_name="lobster").find_one(filter={"basis_id": 1})
        self._check_run(d, mode="lobsternormal")

        wf = self.lp.get_wf_by_fw_id(1)
        self.assertTrue(all([s == 'COMPLETED' for s in wf.fw_states.values()]))

    def test_single_vasp_lobster(self):
        # will only test lobster_calculation_1
        self._single_vasp_lobster(fake=True)

    def test_single_lobster_db_insertion(self):
        self._single_lobster_db_insertion(fake=True)

    # def test_single_vasp_lobster_integration(self):
    #     # integration test
    #     if VASP_CMD and LOBSTER_CMD:
    #         self._single_vasp_lobster(fake=False)


if __name__ == "__main__":
    unittest.main()
