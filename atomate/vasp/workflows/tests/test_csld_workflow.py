# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
from monty.os.path import which
import unittest
import numpy as np
import shutil

from pymatgen import Structure
from atomate.vasp.workflows.base.csld import CompressedSensingLatticeDynamicsWF
from atomate.vasp.firetasks.parse_outputs import (
    CSLDForceConstantsToDB,
    ShengBTEToDB,
)
from atomate.utils.testing import AtomateTest, DB_DIR

import json

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
ref_dir = os.path.join(module_dir, "..", "..", "test_files", "csld_wf")

__author__ = "Rees Chang"
__email__ = "rc564@cornell.edu"

csld_present = which("csld_main")

class TestCompressedSensingLatticeDynamicsWorkflow(AtomateTest):

    @staticmethod
    def init():
        with open(os.path.join(ref_dir, "FW.json"),
                  "r") as f:
            sample_task = json.load(f)[0]
        wf_uuid = sample_task["spec"]["_tasks"][0]["wf_uuid"]
        parent_structure = Structure.from_dict(
            sample_task["spec"]["_tasks"][0]["parent_structure"])
        trans_mat = sample_task["spec"]["_tasks"][0]["trans_mat"]
        supercell_structure = Structure.from_dict(
            sample_task["spec"]["_tasks"][0]["supercell_structure"])
        supercell_smallest_dim = sample_task["spec"]["_tasks"][0]["supercell_smallest_dim"]
        perturbed_supercells = sample_task["spec"]["_tasks"][0]["perturbed_supercells"]
        for i in range(len(perturbed_supercells)):
            perturbed_supercells[i] = Structure.from_dict(
                perturbed_supercells[i])
        disps = sample_task["spec"]["_tasks"][0]["disps"]
        first_pass = sample_task["spec"]["_tasks"][0]["first_pass"]
        static_user_incar_settings = sample_task["spec"]["_tasks"][0]["static_user_incar_settings"]
        env_vars = sample_task["spec"]["_tasks"][0]["env_vars"]
        shengbte_t_range = sample_task["spec"]["_tasks"][0]["shengbte_t_range"]
        shengbte_fworker = sample_task["spec"]["_tasks"][0]["shengbte_fworker"]

        csld_firetask = CSLDForceConstantsToDB(
            db_file=os.path.join(DB_DIR, "db.json"),
            wf_uuid=wf_uuid,
            parent_structure=parent_structure,
            trans_mat=trans_mat,
            supercell_structure=supercell_structure,
            supercell_smallest_dim=supercell_smallest_dim,
            perturbed_supercells=perturbed_supercells,
            disps=disps,
            first_pass=first_pass,
            static_user_incar_settings=static_user_incar_settings,
            env_vars=env_vars,
            shengbte_t_range=shengbte_t_range,
            shengbte_fworker=shengbte_fworker,
            force_diagonal_transformation=True,
        )

        return csld_firetask

    def test_collect_successful_static_calcs_results(self):
        tasks = self.get_task_collection()
        csld_firetask = self.init()
        with open(os.path.join(ref_dir, "FW.json"),
                  "r") as f:
            sample_task = json.load(f)[1]
        tasks.insert_one(sample_task)

        supercells_forces, successful_disps = csld_firetask.collect_successful_static_calcs_results(
            tasks)

        self.assertEqual(
            float(successful_disps[0]),
            0.1
        )
        self.assertTrue(
            isinstance(supercells_forces, list)
        )
        self.assertEqual(
            supercells_forces[0][0][0],
            0.07831190000000000373
        )
        self.assertEqual(
            supercells_forces[0][0][1],
            0.13348752999999999314
        )
        self.assertEqual(
            supercells_forces[0][0][2],
            -0.7885310199999999714
        )
        self.assertEqual(
            len(supercells_forces[0]),
            216
        )

    def test_set_params(self):

        csld_firetask = self.init()
        csld_firetask.set_params(
            iteration_number=0,
            disps=[0.1],
            cluster_diam='11 6.5 5.0',
            max_order=3,
            submodel1='anh 0 1 2 3',
            export_sbte='5 5 5 2 3',
            supercells_forces=[np.eye(3)]
        )

        self.assertEqual(
            csld_firetask["csld_settings"]["structure"]["prim"],
            "Si_supercell_iter0/POSCAR")
        self.assertEqual(
            csld_firetask["csld_settings"]["model"]["max_order"],
            '3'
        )
        self.assertEqual(
            csld_firetask["csld_settings"]["model"]["cluster_diameter"],
            '11 6.5 5.0'
        )
        self.assertEqual(
            csld_firetask["csld_settings"]["model"]["cluster_filter"],
            r"lambda cls: ((cls.order_uniq <= 2) or (cls.bond_counts(2.9) "
            r">= 2)) and cls.is_small_in('" + "Si_supercell_iter0" + "/sc.txt')"
        )
        self.assertEqual(
            csld_firetask["csld_settings"]["training"]["traindat1"],
            'Si_supercell_iter0/SPOSCAR Si_supercell_iter0/disp0.10'
        )
        self.assertEqual(
            csld_firetask["csld_settings"]["fitting"]["submodel1"],
            'anh 0 1 2 3'
        )
        self.assertEqual(
            csld_firetask["csld_settings"]["export_potential"]["export_shengbte"],
            '5 5 5 2 3'
        )
        self.assertEqual(
            csld_firetask["forces_paths"],
            ['Si_supercell_iter0/disp0.10']
        )
        self.assertTrue(
            os.path.isfile(csld_firetask["forces_paths"][0] + "/force.txt")
        )
        self.assertTrue(
            os.path.isfile('Si_supercell_iter0/POSCAR')
        )
        self.assertTrue(
            os.path.isfile('Si_supercell_iter0/SPOSCAR')
        )
        self.assertTrue(
            os.path.isfile('Si_supercell_iter0/sc.txt')
        )
        shutil.rmtree('Si_supercell_iter0')

    @unittest.skipIf(not csld_present, "CSLD not present")
    def test_run_csld(self):
        print('hi')





