# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

from monty.os.path import which

from atomate.vasp.workflows.base.magnetism import MagneticOrderingsWF
from atomate.vasp.firetasks.parse_outputs import (
    MagneticDeformationToDB,
    MagneticOrderingsToDB,
)
from atomate.utils.testing import AtomateTest, DB_DIR

from json import load
from pymatgen import Structure

__author__ = "Matthew Horton"
__email__ = "mkhorton@lbl.gov"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
ref_dir = os.path.join(module_dir, "..", "..", "test_files", "magnetism_wf")

enum_cmd = which("enum.x") or which("multienum.x")
makestr_cmd = which("makestr.x") or which("makeStr.x") or which("makeStr.py")
enumlib_present = enum_cmd and makestr_cmd


class TestMagneticOrderingsWorkflow(AtomateTest):

    @unittest.skipIf(not enumlib_present, "enumlib not present")
    def test_ordering_enumeration(self):

        # simple afm
        structure = Structure.from_file(os.path.join(ref_dir, "ordering/LaMnO3.json"))
        wf = MagneticOrderingsWF(structure)
        self.assertEqual(wf.input_origin, "afm")

        # ferrimagnetic (Cr produces net spin)
        structure = Structure.from_file(os.path.join(ref_dir, "ordering/Cr2NiO4.json"))
        wf = MagneticOrderingsWF(structure)
        self.assertEqual(wf.input_origin, "ferri_by_Cr")

        # antiferromagnetic on single magnetic site
        structure = Structure.from_file(os.path.join(ref_dir, "ordering/Cr2WO6.json"))
        wf = MagneticOrderingsWF(structure)
        self.assertEqual(wf.input_origin, "afm_by_Cr")

        # afm requiring large cell size
        # (enable for further development of workflow, too slow for CI)

        # structure = Structure.from_file(os.path.join(ref_dir, "CuO.json"))
        # wf = MagneticOrderingsWF(structure, default_magmoms={'Cu': 1.73},
        #                         transformation_kwargs={'max_cell_size': 4})
        # self.assertEqual(wf.input_origin, "afm")

        # antiferromagnetic by structural motif
        structure = Structure.from_file(os.path.join(ref_dir, "ordering/Ca3Co2O6.json"))
        wf = MagneticOrderingsWF(
            structure,
            strategies=("antiferromagnetic_by_motif",),
            # this example just misses default cut-off, so do not truncate
            truncate_by_symmetry=False,
            transformation_kwargs={"max_cell_size": 2},
        )
        self.assertEqual(wf.input_origin, "afm_by_motif_2a")

    def test_analysis(self):

        # load example tasks (since workflow re-uses existing FW building
        # blocks for the actual calculations, the most important test is
        # the new analysis task)
        tasks = self.get_task_collection()
        with open(os.path.join(ref_dir, "ordering/sample_tasks.json"), "r") as f:
            sample_tasks = load(f)
        wf_uuid = sample_tasks[0]["wf_meta"]["wf_uuid"]
        parent_structure = Structure.from_dict(sample_tasks[0]["input"]["structure"]).get_primitive_structure()
        tasks.insert_many(sample_tasks)

        toDb = MagneticOrderingsToDB(
            db_file=os.path.join(DB_DIR, "db.json"), wf_uuid=wf_uuid,
            parent_structure=parent_structure,
            perform_bader=False, scan=False
        )
        toDb.run_task({})

        mag_ordering_collection = self.get_task_database().magnetic_orderings
        from pprint import pprint
        stable_ordering = mag_ordering_collection.find_one({"stable": True})
        self.assertEqual(stable_ordering['input']['index'], 2)
        self.assertAlmostEqual(stable_ordering['magmoms']['vasp'][0], -2.738)


class TestMagneticDeformationWorkflow(AtomateTest):
    def test_analysis(self):

        # load example tasks (since workflow re-uses existing FW building
        # blocks for the actual calculations, the most important test is
        # the new analysis task)
        tasks = self.get_task_collection()
        with open(os.path.join(ref_dir, "deformation/sample_tasks.json"), "r") as f:
            sample_tasks = load(f)
        wf_uuid = sample_tasks[0]["wf_meta"]["wf_uuid"]
        tasks.insert_many(sample_tasks)

        toDb = MagneticDeformationToDB(
            db_file=os.path.join(DB_DIR, "db.json"), wf_uuid=wf_uuid
        )
        toDb.run_task({})

        mag_def_collection = self.get_task_database().magnetic_deformation
        magnetic_deformation = mag_def_collection.find_one({"success": True})[
            "magnetic_deformation"
        ]

        self.assertAlmostEqual(magnetic_deformation, 0.12076495213674536)
