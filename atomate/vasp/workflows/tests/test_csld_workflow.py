# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

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


class TestCompressedSensingLatticeDynamicsWorkflow(AtomateTest):

    # def test_CSLD(self):
    #
    #     structure = Structure.from_file(os.path.join(ref_dir, "POSCAR-well_relaxed_Si"))
    #     wf = CompressedSensingLatticeDynamicsWF(structure,
    #                                             symmetrize=False,
    #                                             num_nn_dists=6,
    #                                             num_displacements=10,
    #                                             supercells_per_displacement_distance=1,
    #                                             force_diagonal_transformation=True
    #                                             )
    #
    #     tasks = self.get_task_collection()
    #     with

    def test_set_params(self):
        # tasks = self.get_task_collection()
        with open(os.path.join(ref_dir, "FW.json"),
                  "r") as f:
            sample_task = json.load(f)
        wf_uuid = sample_task["spec"]["_tasks"][0]["wf_uuid"]
        parent_structure = Structure.from_dict(
            sample_task["spec"]["_tasks"][0]["parent_structure"])
        trans_mat = sample_task["spec"]["_tasks"][0]["trans_mat"]
        supercell_structure = sample_task["spec"]["_tasks"][0]["supercell_structure"]
        supercell_smallest_dim = sample_task["spec"]["_tasks"][0]["supercell_smallest_dim"]
        perturbed_supercells = sample_task["spec"]["_tasks"][0]["perturbed_supercells"]
        disps = sample_task["spec"]["_tasks"][0]["disps"]
        first_pass = sample_task["spec"]["_tasks"][0]["first_pass"]
        static_user_incar_settings = sample_task["spec"]["_tasks"][0]["static_user_incar_settings"]
        env_vars = sample_task["spec"]["_tasks"][0]["env_vars"]
        shengbte_t_range = sample_task["spec"]["_tasks"][0]["shengbte_t_range"]
        shengbte_fworker = sample_task["spec"]["_tasks"][0]["shengbte_fworker"]

        print(wf_uuid)
        print(parent_structure)
        print(trans_mat)
        print(supercell_structure)
        print(supercell_smallest_dim)
        print(perturbed_supercells)
        print(disps)
        print(first_pass)
        print(static_user_incar_settings)
        print(env_vars)
        print(shengbte_t_range)
        print(shengbte_fworker)

        # csld_firetask = CSLDForceConstantsToDB(
        #     db_file=os.path.join(DB_DIR, "db.json"),  # wot
        #     wf_uuid=wf_uuid,
        #     parent_structure=parent_structure,
        #     trans_mat=trans_mat,
        #     supercell_structure=supercell_structure,
        #     supercell_smallest_dim=supercell_smallest_dim,
        #     perturbed_supercells=perturbed_supercells,
        #     disps=disps,
        #     first_pass=first_pass,
        #     static_user_incar_settings=static_user_incar_settings,
        #     env_vars=env_vars,
        #     shengbte_t_range=shengbte_t_range,
        #     shengbte_fworker=shengbte_fworker
        # )

