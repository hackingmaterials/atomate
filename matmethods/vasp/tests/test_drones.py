# coding: utf-8
# Copyright (c) Materials Virtual Lab.
# Distributed under the terms of the BSD License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import os
import unittest


from matmethods.vasp.drones import MMVaspToDbTaskDrone


module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class MMVaspToDbTaskDroneTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if not os.environ.get("VASP_PSP_DIR"):
            raise unittest.SkipTest(
                'This system is not set up to run VASP jobs. '
                'Please set your VASP_PSP_DIR environment variable.')

        cls.relax = os.path.join(module_dir, "reference_files",
                                 "Si_structure_optimization", "outputs")

    def test_assimimlate(self):
        drone = MMVaspToDbTaskDrone()
        doc = drone.get_task_doc(self.relax)
        # Only the main changes from the vasprun as dict format and currently
        # used schema in pymatgen-db are tested for now.
        self.assertEqual(doc["composition_reduced"], {'Si': 1.0})
        self.assertEqual(doc["composition_unit_cell"], {'Si': 2.0})
        self.assertAlmostEqual(doc["output"]["energy"], -10.84671647)
        self.assertEqual(doc["formula_pretty"], 'Si')
        self.assertEqual(doc["formula_anonymous"], 'A')


if __name__ == "__main__":
    unittest.main()
