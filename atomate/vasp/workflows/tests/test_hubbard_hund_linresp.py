# coding: utf-8

import os
import unittest

from monty.os.path import which

from atomate.vasp.workflows.base.hubbard_hund_linresp import get_wf_hubbard_hund_linresp
from atomate.vasp.firetasks.parse_outputs import HubbardHundLinRespToDb
from atomate.utils.testing import AtomateTest, DB_DIR

from json import load
from pymatgen import Structure

__author__ = "Guy Moore"
__email__ = "gmoore@lbl.gov"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
ref_dir = os.path.join(module_dir, "..", "..", "test_files", "hubbard_hund_linresp_wf")

class TestHubbardHundLinRespWorkflow(AtomateTest):

    def test_analysis(self):

        # load example tasks (since workflow re-uses existing FW building
        # blocks for the actual calculations, the most important test is
        # the new analysis task)

        # Example: LiNiPO4
        formula_pretty = "LiNiPO4"
        tasks = self.get_task_collection()
        with open(os.path.join(ref_dir, "sample_tasks_multisite_spinpolarized.json"), "r") as f:
            sample_tasks = load(f)
        wf_uuid = sample_tasks[0]["wf_meta"]["wf_uuid"]
        parent_structure = Structure.from_dict(sample_tasks[0]["input"]["structure"])
        tasks.insert_many(sample_tasks)

        num_perturb = 2
        spin_polarized = True
        relax_nonmagnetic = True
        toDb = HubbardHundLinRespToDb(
            num_perturb=num_perturb,
            spin_polarized=spin_polarized,
            relax_nonmagnetic=relax_nonmagnetic,
            db_file=os.path.join(DB_DIR, "db.json"), wf_uuid=wf_uuid,
        )
        toDb.run_task({})

        hubbard_hund_linresp_collection = self.get_task_database().hubbard_hund_linresp
        from pprint import pprint
        hubbard_hund_analysis = hubbard_hund_linresp_collection.find_one({"formula_pretty": formula_pretty})

        # For LiNiPO4: site0 = O-p, site1 = Ni-d
        uj_dict = hubbard_hund_analysis["hubbard_hund_results"]

        print("output = ", 
              uj_dict["point"]["values"]["site0"]["U"]["value"],
              uj_dict["point"]["values"]["site1"]["U"]["value"],
              uj_dict["atom"]["values"]["site0"]["simple"]["U"]["value"],
              uj_dict["atom"]["values"]["site1"]["simple"]["U"]["value"],
              uj_dict["atom"]["values"]["site1"]["simple"]["J"]["value"],
              uj_dict["full"]["values"]["site0"]["simple"]["U"]["value"],
              uj_dict["full"]["values"]["site1"]["simple"]["U"]["value"],
              uj_dict["full"]["values"]["site1"]["simple"]["J"]["value"])

        # Point-wise inversion
        self.assertAlmostEqual(uj_dict["point"]["values"]["site0"]["simple"]["U"]["value"],
                               2.278473, 4)
        self.assertAlmostEqual(uj_dict["point"]["values"]["site1"]["simple"]["U"]["value"],
                               2.745100, 4)

        # Atom-wise inversion
        self.assertAlmostEqual(uj_dict["atom"]["values"]["site0"]["simple"]["U"]["value"],
                               9.108129, 4)
        self.assertAlmostEqual(uj_dict["atom"]["values"]["site1"]["simple"]["U"]["value"],
                               4.698635, 4)


        self.assertAlmostEqual(uj_dict["atom"]["values"]["site1"]["simple"]["J"]["value"],
                               0.402767, 4)

        # Full inversion
        self.assertAlmostEqual(uj_dict["full"]["values"]["site0"]["simple"]["U"]["value"],
                               9.096978, 4)
        self.assertAlmostEqual(uj_dict["full"]["values"]["site1"]["simple"]["U"]["value"],
                               4.687404, 4)


        self.assertAlmostEqual(uj_dict["full"]["values"]["site1"]["simple"]["J"]["value"],
                               0.404227, 4)
