# coding: utf-8


import json
import os
import unittest
import zlib

import gridfs
from pymongo import DESCENDING

from fireworks import FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.vasp.powerups import use_custodian, add_namefile, use_fake_vasp, add_trackers, add_bandgap_check, use_potcar_spec
from atomate.vasp.workflows.base.core import get_wf
from atomate.utils.testing import AtomateTest
from atomate.vasp.firetasks.parse_outputs import VaspDrone
from atomate.vasp.database import VaspCalcDb


from pymatgen.io.vasp import Incar
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet, MPScanRelaxSet
from pymatgen.util.testing import PymatgenTest
from pymatgen.core import Structure

__author__ = 'Anubhav Jain, Kiran Mathew'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
reference_dir = os.path.join(module_dir, "..", "..", "test_files")

ref_dirs_si = {"structure optimization": os.path.join(reference_dir, "Si_structure_optimization"),
             "static": os.path.join(reference_dir, "Si_static"),
             "nscf uniform": os.path.join(reference_dir, "Si_nscf_uniform"),
             "nscf line": os.path.join(reference_dir, "Si_nscf_line")}

_fworker = FWorker(env={"db_file": os.path.join(db_dir, "db.json")})

DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...


class TestVaspWorkflows(AtomateTest):

    def setUp(self):
        super(TestVaspWorkflows, self).setUp()
        self.struct_si = PymatgenTest.get_structure("Si")

    def _check_run(self, d, mode):
        if mode not in ["structure optimization", "static", "nscf uniform",
                        "nscf line", "additional field"]:
            raise ValueError("Invalid mode!")

        self.assertEqual(d["formula_pretty"], "Si")
        self.assertEqual(d["formula_anonymous"], "A")
        self.assertEqual(d["nelements"], 1)
        self.assertEqual(d["state"], "successful")
        self.assertAlmostEqual(d["calcs_reversed"][0]["output"]["structure"]["lattice"]["a"], 3.867, 2)
        self.assertEqual(d["output"]["is_gap_direct"], False)

        if mode in ["structure optimization", "static"]:
            self.assertAlmostEqual(d["output"]["energy"], -10.850, 2)
            self.assertAlmostEqual(d["output"]["energy_per_atom"], -5.425, 2)

        if mode == "additional field":
            self.assertAlmostEqual(d["test_additional_field"]["lattice"]["a"], 3.8401979337)

        elif mode in ["ncsf uniform"]:
            self.assertAlmostEqual(d["output"]["energy"], -10.828, 2)
            self.assertAlmostEqual(d["output"]["energy_per_atom"], -5.414, 2)

        self.assertAlmostEqual(d["output"]["bandgap"], 0.65, 1)

        if "nscf" in mode:
            self.assertEqual(d["calcs_reversed"][0]["output"]["outcar"]["total_magnetization"], None)
        else:
            self.assertAlmostEqual(d["calcs_reversed"][0]["output"]["outcar"]["total_magnetization"], 0, 3)

        self.assertLess(d["run_stats"]["overall"]["Elapsed time (sec)"], 180)  # run should take under 3 minutes

        # check the DOS and band structure
        if mode == "nscf uniform" or mode == "nscf line":
            fs = gridfs.GridFS(self.get_task_database(), 'bandstructure_fs')

            # check the band structure
            bs_fs_id = d["calcs_reversed"][0]["bandstructure_fs_id"]
            bs_json = zlib.decompress(fs.get(bs_fs_id).read())
            bs = json.loads(bs_json.decode())
            self.assertEqual(bs["is_spin_polarized"], False)
            self.assertEqual(bs["band_gap"]["direct"], False)
            self.assertAlmostEqual(bs["band_gap"]["energy"], 0.65, 1)
            self.assertEqual(bs["is_metal"], False)

            if mode == "nscf uniform":
                for k in ["is_spin_polarized", "band_gap", "structure",
                          "kpoints", "is_metal", "vbm", "cbm", "labels_dict",
                          "projections", "lattice_rec", "bands"]:
                    self.assertTrue(k in bs)
                    self.assertIsNotNone(bs[k])

                self.assertEqual(bs["@class"], "BandStructure")

            else:
                for k in ["is_spin_polarized", "band_gap", "structure",
                          "kpoints", "is_metal", "vbm", "cbm", "labels_dict",
                          "projections", "lattice_rec", "bands", "branches"]:
                    self.assertTrue(k in bs)
                    self.assertIsNotNone(bs[k])
                self.assertEqual(bs["@class"], "BandStructureSymmLine")

            # check the DOS
            if mode == "nscf uniform":
                fs = gridfs.GridFS(self.get_task_database(), 'dos_fs')
                dos_fs_id = d["calcs_reversed"][0]["dos_fs_id"]

                dos_json = zlib.decompress(fs.get(dos_fs_id).read())
                dos = json.loads(dos_json.decode())
                for k in ["densities", "energies", "pdos", "spd_dos", "atom_dos", "structure"]:
                    self.assertTrue(k in dos)
                    self.assertIsNotNone(dos[k])

                self.assertAlmostEqual(dos["spd_dos"]["p"]["efermi"], 5.625, 1)
                self.assertAlmostEqual(dos["atom_dos"]["Si"]["efermi"], 5.625, 1)
                self.assertAlmostEqual(dos["structure"]["lattice"]["a"], 3.867, 2)
                self.assertAlmostEqual(dos["spd_dos"]["p"]["efermi"], 5.625, 1)
                self.assertAlmostEqual(dos["atom_dos"]["Si"]["efermi"], 5.625, 1)
                self.assertAlmostEqual(dos["structure"]["lattice"]["a"], 3.867, 2)

    def test_single_Vasp(self):
        # add the workflow
        structure = self.struct_si
        my_wf = get_wf(structure, "optimize_only.yaml", vis=MPRelaxSet(structure, force_gamma=True),
                       common_params={"vasp_cmd": VASP_CMD})
        if not VASP_CMD:
            my_wf = use_fake_vasp(my_wf, ref_dirs_si)
        else:
            my_wf = use_custodian(my_wf)
        self.lp.add_wf(my_wf)

        # run the workflow
        rapidfire(self.lp, fworker=_fworker)

        d = self.get_task_collection().find_one({"task_label": "structure optimization"})
        self._check_run(d, mode="structure optimization")

        wf = self.lp.get_wf_by_fw_id(1)
        self.assertTrue(all([s == 'COMPLETED' for s in wf.fw_states.values()]))

    def test_single_Vasp_dbinsertion(self):
        # add the workflow
        structure = self.struct_si
        # instructs to use db_file set by FWorker, see env_chk
        my_wf = get_wf(structure, "optimize_only.yaml", vis=MPRelaxSet(structure, force_gamma=True),
                       common_params={"vasp_cmd": VASP_CMD,
                                      "db_file": ">>db_file<<"})
        if not VASP_CMD:
            my_wf = use_fake_vasp(my_wf, ref_dirs_si)
        else:
            my_wf = use_custodian(my_wf)

        # add an msonable object to additional fields
        my_wf.fws[0].tasks[-1]['additional_fields'].update(
            {"test_additional_field": self.struct_si})
        self.lp.add_wf(my_wf)

        # run the workflow
        rapidfire(self.lp, fworker=_fworker)

        d = self.get_task_collection().find_one()
        self._check_run(d, mode="structure optimization")
        self._check_run(d, mode="additional field")

        wf = self.lp.get_wf_by_fw_id(1)
        self.assertTrue(all([s == 'COMPLETED' for s in wf.fw_states.values()]))

    def test_bandstructure_Vasp(self):
        # add the workflow
        structure = self.struct_si
        # instructs to use db_file set by FWorker, see env_chk
        my_wf = get_wf(structure, "bandstructure.yaml",
                       vis=MPRelaxSet(structure, force_gamma=True),
                       common_params={"vasp_cmd": VASP_CMD,
                                      "db_file": ">>db_file<<"})
        if not VASP_CMD:
            my_wf = use_fake_vasp(my_wf, ref_dirs_si)
        else:
            my_wf = use_custodian(my_wf)

        my_wf = add_namefile(my_wf)  # add a slug of fw-name to output files

        self.lp.add_wf(my_wf)

        # run the workflow
        # set the db_file variable
        rapidfire(self.lp, fworker=_fworker)

        # make sure the structure relaxation ran OK
        d = self.get_task_collection().find_one({"task_label": "structure optimization"},
                                                sort=[("_id", DESCENDING)])
        self._check_run(d, mode="structure optimization")

        # make sure the static run ran OK
        d = self.get_task_collection().find_one({"task_label": "static"}, sort=[("_id", DESCENDING)])
        self._check_run(d, mode="static")

        # make sure the uniform run ran OK
        d = self.get_task_collection().find_one({"task_label": "nscf uniform"}, sort=[("_id", DESCENDING)])
        self._check_run(d, mode="nscf uniform")

        # make sure the uniform run ran OK
        d = self.get_task_collection().find_one({"task_label": "nscf line"}, sort=[("_id", DESCENDING)])
        self._check_run(d, mode="nscf line")

        wf = self.lp.get_wf_by_fw_id(1)
        self.assertTrue(all([s == 'COMPLETED' for s in wf.fw_states.values()]))

    def test_bandgap_check_Vasp(self):
        # add the workflow
        structure = self.struct_si
        # instructs to use db_file set by FWorker, see env_chk
        my_wf = get_wf(structure, "bandstructure.yaml",
                       vis=MPRelaxSet(structure, force_gamma=True),
                       common_params={"vasp_cmd": VASP_CMD,
                                      "db_file": ">>db_file<<"})
        if not VASP_CMD:
            my_wf = use_fake_vasp(my_wf, ref_dirs_si)
        else:
            my_wf = use_custodian(my_wf)

        my_wf = add_namefile(my_wf)  # add a slug of fw-name to output files
        my_wf = add_bandgap_check(my_wf, check_bandgap_params={"max_gap": 0.1}, fw_name_constraint="structure optimization")
        self.lp.add_wf(my_wf)

        # run the workflow
        # set the db_file variable
        rapidfire(self.lp, fworker=_fworker)

        # structure optimization should be completed
        self.assertEqual(self.lp.fireworks.find_one(
            {"name": "Si-structure optimization"}, {"state": 1})["state"],
                         "COMPLETED")

        self.assertEqual(self.lp.fireworks.find_one(
            {"name": "Si-static"}, {"state": 1})["state"],
                         "DEFUSED")

    def test_trackers(self):
        # add the workflow
        structure = self.struct_si
        my_wf = get_wf(structure, "optimize_only.yaml",
                       vis=MPRelaxSet(structure, force_gamma=True),
                       common_params={"vasp_cmd": VASP_CMD})

        if not VASP_CMD:
            my_wf = use_fake_vasp(my_wf, ref_dirs_si)
        else:
            my_wf = use_custodian(my_wf)

        my_wf = add_trackers(my_wf)
        self.lp.add_wf(my_wf)

        # run the workflow
        rapidfire(self.lp, fworker=_fworker)

        for x in self.lp.get_tracker_data(1):
            for t in x["trackers"]:
                self.assertGreater(len(t.content.split("\n")), 20)

        wf = self.lp.get_wf_by_fw_id(1)
        self.assertTrue(all([s == 'COMPLETED' for s in wf.fw_states.values()]))


    def test_chgcar_db_read_write(self):
        # generate a doc from the test folder
        drone = VaspDrone(parse_chgcar=True, parse_aeccar=True)
        print(ref_dirs_si['static'])
        doc = drone.assimilate(ref_dirs_si['static']+'/outputs')
        # insert the doc make sure that the
        cc = doc['calcs_reversed'][0]['chgcar']
        self.assertAlmostEqual(cc.data['total'].sum()/cc.ngridpts, 8.0, 4)
        cc = doc['calcs_reversed'][0]['aeccar0']
        self.assertAlmostEqual(cc.data['total'].sum()/cc.ngridpts, 23.253588293583313, 4)
        cc = doc['calcs_reversed'][0]['aeccar2']
        self.assertAlmostEqual(cc.data['total'].sum()/cc.ngridpts, 8.01314480789829, 4)
        mmdb = VaspCalcDb.from_db_file(os.path.join(db_dir, "db.json"))
        t_id = mmdb.insert_task(doc, use_gridfs=True)
        # space is freed up after uploading the document
        self.assertRaises(KeyError, lambda: doc['calcs_reversed'][0]['chgcar'])
        self.assertRaises(KeyError, lambda: doc['calcs_reversed'][0]['aeccar0'])
        self.assertRaises(KeyError, lambda: doc['calcs_reversed'][0]['aeccar2'])
        cc = mmdb.get_chgcar(task_id=t_id)
        self.assertAlmostEqual(cc.data['total'].sum()/cc.ngridpts, 8.0, 4)
        dcc = mmdb.get_aeccar(task_id=t_id)
        self.assertAlmostEqual(dcc['aeccar0'].data['total'].sum()/cc.ngridpts, 23.253588293583313, 4)
        self.assertAlmostEqual(dcc['aeccar2'].data['total'].sum()/cc.ngridpts, 8.01314480789829, 4)
        # check the retrieve_task function for the same fake calculation
        ret_task = mmdb.retrieve_task(t_id)
        ret_chgcar = ret_task['calcs_reversed'][0]['chgcar']
        ret_aeccar0 = ret_task['calcs_reversed'][0]['aeccar0']
        ret_aeccar2 = ret_task['calcs_reversed'][0]['aeccar2']
        ret_aeccar = ret_aeccar0 + ret_aeccar2
        self.assertAlmostEqual(ret_chgcar.data['total'].sum()/ret_chgcar.ngridpts, 8.0, 4)
        self.assertAlmostEqual(ret_aeccar.data['total'].sum()/ret_aeccar.ngridpts, 31.2667331015, 4)

    def test_chgcar_db_read(self):
        # add the workflow
        structure = self.struct_si
        # instructs to use db_file set by FWorker, see env_chk
        my_wf = get_wf(structure, "static_only.yaml", vis=MPStaticSet(structure, force_gamma=True),
                       common_params={"vasp_cmd": VASP_CMD,
                                      "db_file": ">>db_file<<"})
        if not VASP_CMD:
            my_wf = use_fake_vasp(my_wf, ref_dirs_si)
        else:
            my_wf = use_custodian(my_wf)

        # set the flags for storing charge densties
        my_wf.fws[0].tasks[-1]["parse_chgcar"] = True
        my_wf.fws[0].tasks[-1]["parse_aeccar"] = True
        self.lp.add_wf(my_wf)

        # run the workflow
        # set the db_file variable
        rapidfire(self.lp, fworker=_fworker)

        d = self.get_task_collection().find_one()
        self._check_run(d, mode="static")

        wf = self.lp.get_wf_by_fw_id(1)
        self.assertTrue(all([s == 'COMPLETED' for s in wf.fw_states.values()]))

        chgcar_fs_id = d["calcs_reversed"][0]["chgcar_fs_id"]
        accar0_fs_id = d["calcs_reversed"][0]["aeccar0_fs_id"]
        accar2_fs_id = d["calcs_reversed"][0]["aeccar2_fs_id"]

        self.assertTrue(bool(chgcar_fs_id))
        self.assertTrue(bool(accar0_fs_id))
        self.assertTrue(bool(accar2_fs_id))


class TestScanOptimizeWorkflow(AtomateTest):

    def setUp(self):
        super(TestScanOptimizeWorkflow, self).setUp()

    def _run_scan_relax(self, wf, dir_name):
        if not VASP_CMD:
            wf = use_fake_vasp(wf,
                               {"SCAN structure optimization": os.path.join(
                                reference_dir, dir_name)},
                               check_kpoints=False,
                               check_potcar=False,
                               clear_inputs=False,
                               check_incar=False
                               )
        else:
            wf = use_custodian(wf)

        wf = use_potcar_spec(wf)
        self.lp.add_wf(wf)

        # run the workflow
        rapidfire(self.lp, fworker=_fworker)

    def _get_launch_dir(self):
        # retrieve the launcher directory
        d = list(self.get_task_collection().find({"task_label": "SCAN structure optimization"}))[-1]
        launch_dir = d["dir_name"].split(":")[1]
        return launch_dir

    def test_SCAN_no_bandgap(self):
        # A structure with bandgap = 0 (default) should have KSPACING equal to 0.22
        structure = Structure.from_file(os.path.join(reference_dir, "SCAN_structure_optimization_Al/inputs", "POSCAR"))

        my_wf = get_wf(structure, "SCAN_optimization.yaml", vis=MPScanRelaxSet(structure),
                       common_params={"vasp_cmd": VASP_CMD})

        self._run_scan_relax(my_wf, "SCAN_structure_optimization_Al")

        # Check INCAR.orig
        incar_orig = Incar.from_file(os.path.join(self._get_launch_dir(), "INCAR.orig.gz"))
        ref_incar = Incar.from_file(os.path.join(reference_dir, "SCAN_structure_optimization_Al/inputs", "INCAR.orig"))
        for p in incar_orig.keys():
            if p == "MAGMOM":  # Ignore MAGMOM b/c structure initialized from POSCAR cannot have a MAGMOM
                pass
            else:
                self.assertEqual(incar_orig[p], ref_incar[p])

        # Check INCAR.relax1
        incar = Incar.from_file(os.path.join(self._get_launch_dir(), "INCAR.relax1.gz"))
        self.assertEqual(incar["METAGGA"], "None")
        self.assertEqual(incar["EDIFFG"], -0.05)
        self.assertEqual(incar["LWAVE"], False)

        # Check INCAR.relax2
        incar = Incar.from_file(os.path.join(self._get_launch_dir(), "INCAR.relax2.gz"))
        self.assertEqual(incar["METAGGA"], "None")
        self.assertEqual(incar["LWAVE"], True)
        self.assertEqual(incar["NSW"], 0)
        self.assertEqual(incar["EDIFFG"], -0.05)
        self.assertEqual(incar["ICHARG"], 1)
        self.assertEqual(incar["ISTART"], 0)

        # Check INCAR.relax3
        incar = Incar.from_file(os.path.join(self._get_launch_dir(), "INCAR.relax3.gz"))
        for p in incar.keys():
            if p == "KSPACING":
                self.assertEqual(incar[p], 0.22)
            elif p == "ICHARG" or p == "ISTART":
                self.assertEqual(incar[p], 1)
            else:
                self.assertEqual(incar_orig[p], incar[p])

    def test_SCAN_small_bandgap(self):
        # A structure with a small bandgap (LiH) should result in a KSPACING
        # value of 0.351275

        structure = Structure.from_file(os.path.join(reference_dir, "SCAN_structure_optimization_LiH/inputs", "POSCAR"))

        my_wf = get_wf(structure, "SCAN_optimization.yaml", vis=MPScanRelaxSet(structure),
                       common_params={"vasp_cmd": VASP_CMD})

        self._run_scan_relax(my_wf, "SCAN_structure_optimization_LiH")

        # Check INCAR.orig
        incar_orig = Incar.from_file(os.path.join(self._get_launch_dir(), "INCAR.orig.gz"))
        ref_incar = Incar.from_file(os.path.join(reference_dir, "SCAN_structure_optimization_LiH/inputs", "INCAR.orig"))
        for p in incar_orig.keys():
            if p == "MAGMOM":  # Ignore MAGMOM b/c structure initialized from POSCAR cannot have a MAGMOM
                pass
            else:
                self.assertEqual(incar_orig[p], ref_incar[p])

        # Check INCAR.relax1
        incar = Incar.from_file(os.path.join(self._get_launch_dir(), "INCAR.relax1.gz"))
        self.assertEqual(incar["METAGGA"], "None")
        self.assertEqual(incar["EDIFFG"], -0.05)
        self.assertEqual(incar["LWAVE"], False)

        # Check INCAR.relax2
        incar = Incar.from_file(os.path.join(self._get_launch_dir(), "INCAR.relax2.gz"))
        self.assertEqual(incar["METAGGA"], "None")
        self.assertEqual(incar["LWAVE"], True)
        self.assertEqual(incar["NSW"], 0)
        self.assertEqual(incar["EDIFFG"], -0.05)
        self.assertEqual(incar["ICHARG"], 1)
        self.assertEqual(incar["ISTART"], 0)

        # Check INCAR.relax3
        incar = Incar.from_file(os.path.join(self._get_launch_dir(), "INCAR.relax3.gz"))
        for p in incar.keys():
            if p == "KSPACING":
                self.assertAlmostEqual(incar[p], 0.351275, 4)
            elif p == "ICHARG" or p == "ISTART":
                self.assertEqual(incar[p], 1)
            elif p == "ISMEAR":
                self.assertEqual(incar[p], -5)
            elif p == "SIGMA":
                self.assertEqual(incar[p], 0.05)
            else:
                self.assertEqual(incar_orig[p], incar[p])

    def test_SCAN_large_bandgap(self):
        # A structure with a large bandgap (LiF) should result in KSPACING
        # hitting the maximum allowed value of 0.44

        structure = Structure.from_file(os.path.join(reference_dir, "SCAN_structure_optimization_LiF/inputs", "POSCAR"))

        my_wf = get_wf(structure, "SCAN_optimization.yaml", vis=MPScanRelaxSet(structure),
                       common_params={"vasp_cmd": VASP_CMD})

        self._run_scan_relax(my_wf, "SCAN_structure_optimization_LiF")

        # Check INCAR.orig generated by the InputSet
        incar_orig = Incar.from_file(os.path.join(self._get_launch_dir(), "INCAR.orig.gz"))
        ref_incar = Incar.from_file(os.path.join(reference_dir, "SCAN_structure_optimization_LiF/inputs", "INCAR.orig"))
        for p in incar_orig.keys():
            if p == "MAGMOM":  # Ignore MAGMOM b/c structure initialized from POSCAR cannot have a MAGMOM
                pass
            else:
                self.assertEqual(incar_orig[p], ref_incar[p])

        # Check INCAR.relax1 generated by the Workflow
        incar = Incar.from_file(os.path.join(self._get_launch_dir(), "INCAR.relax1.gz"))
        self.assertEqual(incar["METAGGA"], "None")
        self.assertEqual(incar["EDIFFG"], -0.05)
        self.assertEqual(incar["LWAVE"], False)

        # Check INCAR.relax2
        incar = Incar.from_file(os.path.join(self._get_launch_dir(), "INCAR.relax2.gz"))
        self.assertEqual(incar["METAGGA"], "None")
        self.assertEqual(incar["LWAVE"], True)
        self.assertEqual(incar["NSW"], 0)
        self.assertEqual(incar["EDIFFG"], -0.05)
        self.assertEqual(incar["ICHARG"], 1)
        self.assertEqual(incar["ISTART"], 0)


        # Check INCAR.relax3 for the correct kspacing
        incar = Incar.from_file(os.path.join(self._get_launch_dir(), "INCAR.relax3.gz"))
        for p in incar.keys():
            if p == "KSPACING":
                self.assertEqual(incar[p], 0.44)
            elif p == "ICHARG" or p == "ISTART":
                self.assertEqual(incar[p], 1)
            elif p == "ISMEAR":
                self.assertEqual(incar[p], -5)
            elif p == "SIGMA":
                self.assertEqual(incar[p], 0.05)
            else:
                self.assertEqual(incar_orig[p], incar[p])

    def test_SCAN_with_vdw(self):
        # Verify appropriate changes to the INCAR when VdW is enabled
        # VdW should be off for relax1 (GGA) and re-enabled for relax2 (SCAN)

        structure = Structure.from_file(os.path.join(reference_dir, "SCAN_structure_optimization_LiF_vdw/inputs", "POSCAR"))

        my_wf = get_wf(structure, "SCAN_optimization.yaml", vis=MPScanRelaxSet(structure, vdw="rvv10"),
                       common_params={"vasp_cmd": VASP_CMD, "vdw_kernel_dir": os.path.join(reference_dir, "SCAN_structure_optimization_LiF_vdw/inputs")})

        self._run_scan_relax(my_wf, "SCAN_structure_optimization_LiF_vdw")

        # Check INCAR.orig
        incar_orig = Incar.from_file(os.path.join(self._get_launch_dir(), "INCAR.orig.gz"))
        ref_incar = Incar.from_file(os.path.join(reference_dir, "SCAN_structure_optimization_LiF_vdw/inputs", "INCAR.orig"))
        for p in incar_orig.keys():
            if p == "MAGMOM":  # Ignore MAGMOM b/c structure initialized from POSCAR cannot have a MAGMOM
                pass
            else:
                self.assertEqual(incar_orig[p], ref_incar[p])

        # Check INCAR.relax1
        incar = Incar.from_file(os.path.join(self._get_launch_dir(), "INCAR.relax1.gz"))
        self.assertIsNone(incar.get("LUSE_VDW", None))
        self.assertIsNone(incar.get("BPARAM", None))
        self.assertEqual(incar["METAGGA"], "None")
        self.assertEqual(incar["EDIFFG"], -0.05)
        self.assertEqual(incar["LWAVE"], False)

        # Check INCAR.relax2
        incar = Incar.from_file(os.path.join(self._get_launch_dir(), "INCAR.relax2.gz"))
        self.assertEqual(incar["METAGGA"], "None")
        self.assertEqual(incar["LWAVE"], True)
        self.assertEqual(incar["NSW"], 0)
        self.assertEqual(incar["EDIFFG"], -0.05)
        self.assertEqual(incar["ICHARG"], 1)
        self.assertEqual(incar["ISTART"], 0)

        # Check INCAR.relax3
        incar = Incar.from_file(os.path.join(self._get_launch_dir(), "INCAR.relax3.gz"))
        for p in incar.keys():
            if p == "KSPACING":
                self.assertEqual(incar[p], 0.44)
            elif p == "ICHARG" or p == "ISTART":
                self.assertEqual(incar[p], 1)
            elif p == "ISMEAR":
                self.assertEqual(incar[p], -5)
            elif p == "SIGMA":
                self.assertEqual(incar[p], 0.05)
            elif p == "MAGMOM":  # Ignore MAGMOM b/c structure initialized from POSCAR cannot have a MAGMOM
                pass
            else:
                self.assertEqual(incar_orig[p], incar[p])

    def test_SCAN_incar_override(self):
        # user incar settings should be passed all the way through the workflow

        structure = Structure.from_file(os.path.join(reference_dir, "SCAN_structure_optimization_LiH/inputs", "POSCAR"))

        my_wf = get_wf(structure, "SCAN_optimization.yaml",
                       vis=MPScanRelaxSet(structure,
                                          user_potcar_functional="PBE_52",
                                          user_incar_settings={"NSW": 10, "SYMPREC": 1e-6, "SIGMA": 0.1}
                                          ),
                       common_params={"vasp_cmd": VASP_CMD})

        self._run_scan_relax(my_wf, "SCAN_structure_optimization_LiH")

        # Check INCAR.orig
        incar = Incar.from_file(os.path.join(self._get_launch_dir(), "INCAR.orig.gz"))
        self.assertEqual(incar["NSW"], 10)
        self.assertEqual(incar["SYMPREC"], 1e-6)
        self.assertEqual(incar["SIGMA"], 0.1)

        # Check INCAR.relax1
        incar = Incar.from_file(os.path.join(self._get_launch_dir(), "INCAR.relax1.gz"))
        self.assertEqual(incar["NSW"], 10)
        self.assertEqual(incar["SYMPREC"], 1e-6)
        self.assertEqual(incar["SIGMA"], 0.1)

        # Check INCAR.relax2
        incar = Incar.from_file(os.path.join(self._get_launch_dir(), "INCAR.relax2.gz"))
        self.assertEqual(incar["NSW"], 0)
        self.assertEqual(incar["SYMPREC"], 1e-6)
        self.assertEqual(incar["SIGMA"], 0.1)

        # Check INCAR.relax3
        incar = Incar.from_file(os.path.join(self._get_launch_dir(), "INCAR.relax3.gz"))
        self.assertEqual(incar["NSW"], 10)
        self.assertEqual(incar["SYMPREC"], 1e-6)
        self.assertEqual(incar["SIGMA"], 0.1)


if __name__ == "__main__":
    unittest.main()
