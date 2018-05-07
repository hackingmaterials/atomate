# coding: utf-8
# Copyright (c) Materials Virtual Lab.
# Distributed under the terms of the BSD License.

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
from atomate.qchem.drones import QChemDrone
import numpy as np

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class QChemDroneTest(unittest.TestCase):

    # def test_assimilate(self):
    #     drone = QChemDrone()
    #     doc = drone.assimilate(path=os.path.join(module_dir, "..", "test_files", "FF_working"), input_file="test.qin_opt1", output_file="test.qout_opt1")
    #     self.assertEqual(doc[],)

    def test_runs_assimilate(self):
        drone = QChemDrone(runs=['opt_0','freq_0','opt_1','freq_1','opt_2','freq_2','opt_3','freq_3'], additional_fields={"special_run_type": "frequency_flattener"})
        doc = drone.assimilate(path=os.path.join(module_dir, "..", "test_files", "FF_working"), input_file="test.qin", output_file="test.qout")
        self.assertEqual(doc["special_run_type"],"frequency_flattener")
        self.assertEqual(doc["schema"],{'code': 'atomate', 'version': '0.7.5'})
        self.assertEqual(doc["input"]["job_type"],"opt")
        self.assertEqual(doc["output"]["job_type"],"freq")
        self.assertEqual(doc["output"]["frequencies"],array([  12.52,   45.28,  260.96,  329.08,  531.01,  582.11,  744.91, 779.2 ,  800.47,  863.15,  928.68,  969.  , 1092.86, 1124.  ,1147.64, 1209.1 , 1387.39, 1693.97, 1913.05, 3316.2 , 3341.73]))
        self.assertEqual(doc["num_frequencies_flattened"],1.0)
        self.assertEqual(doc["walltime"],935.29)
        self.assertEqual(doc["cputime"],3616.6400000000003)
        self.assertEqual(doc["smiles"],"O1[C](O[Li])OC=C1")
        self.assertIn("custodian", doc)
        self.assertIn("calcs_reversed", doc)
        self.assertIn("initial_molecule", doc["input"])
        self.assertIn("initial_molecule", doc["output"])
        self.assertIn("optimized_molecule", doc["output"])
        self.assertIn("last_updated", doc)
        self.assertIn("dir_name", doc)
        self.assertEqual(len(doc["calcs_reversed"]),4)
        self.assertEqual(doc["output"]["frequencies"],doc["calcs_reversed"][0]["frequencies"])
        self.assertEqual(list(doc["calcs_reversed"][0].keys()),
                         list(doc["calcs_reversed"][2].keys()))
        self.assertEqual(list(doc["calcs_reversed"][1].keys()),
                         list(doc["calcs_reversed"][3].keys()))
        
        # self.assertIn("bandstructure", doc["calcs_reversed"][0])
        # oszicar2 = Oszicar(os.path.join(self.relax2, "OSZICAR.relax2.gz"))
        # outcar1 = Outcar(os.path.join(self.relax2, "OUTCAR.relax1.gz"))
        # outcar2 = Outcar(os.path.join(self.relax2, "OUTCAR.relax2.gz"))
        # outcar1 = outcar1.as_dict()
        # outcar2 = outcar2.as_dict()
        # run_stats1 = outcar1.pop("run_stats")
        # run_stats2 = outcar2.pop("run_stats")
        # self.assertEqual(len(doc["calcs_reversed"]), 2)
        # self.assertEqual(doc["composition_reduced"], {'Si': 1.0})
        # self.assertEqual(doc["composition_unit_cell"], {'Si': 2.0})
        # self.assertAlmostEqual(doc["output"]["energy"], oszicar2.ionic_steps[-1]["E0"])
        # self.assertEqual(doc["formula_pretty"], 'Si')
        # self.assertEqual(doc["formula_anonymous"], 'A')
        # self.assertEqual(list(doc["calcs_reversed"][0]["input"].keys()),
        #                  list(doc["calcs_reversed"][1]["input"].keys()))
        # self.assertEqual(list(doc["calcs_reversed"][0]["output"].keys()),
        #                  list(doc["calcs_reversed"][1]["output"].keys()))
        # self.assertEqual(doc["calcs_reversed"][0]["output"]["energy"], doc["output"]["energy"])
        # self.assertEqual(doc["run_stats"][doc["calcs_reversed"][0]["task"]["name"]], run_stats2)
        # self.assertEqual(doc["run_stats"][doc["calcs_reversed"][1]["task"]["name"]], run_stats1)
        # self.assertEqual(doc["calcs_reversed"][0]["output"]["outcar"], outcar2)
        # self.assertEqual(doc["calcs_reversed"][1]["output"]["outcar"], outcar1)


if __name__ == "__main__":
    unittest.main()
