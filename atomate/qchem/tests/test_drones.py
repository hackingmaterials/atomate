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

    def test_assimilate_opt(self):
        drone = QChemDrone()
        doc = drone.assimilate(path=os.path.join(module_dir, "..", "test_files", "FF_working"), input_file="test.qin_opt1", output_file="test.qout_opt1")
        self.assertEqual(doc["input"]["job_type"],"opt")
        self.assertEqual(doc["output"]["job_type"],"opt")
        self.assertEqual(doc["output"]["final_energy"],-348.652462579636)
        self.assertEqual(doc["walltime"],62.83)
        self.assertEqual(doc["cputime"],715.76)
        self.assertEqual(doc["smiles"],"O1[C](O[Li])OC=C1")
        self.assertEqual(doc["formula_pretty"],"LiH2(CO)3")
        self.assertEqual(doc["formula_anonymous"],"AB2C3D3")
        self.assertEqual(doc["chemsys"],"C-H-Li-O")
        self.assertEqual(doc["pointgroup"],"Cs")
        self.assertIn("custodian", doc)
        self.assertIn("calcs_reversed", doc)
        self.assertIn("initial_molecule", doc["input"])
        self.assertIn("initial_molecule", doc["output"])
        self.assertIn("optimized_molecule", doc["output"])
        self.assertIn("last_updated", doc)
        self.assertIn("dir_name", doc)
        self.assertEqual(len(doc["calcs_reversed"]),1)

    def test_assimilate_freq(self):
        drone = QChemDrone()
        doc = drone.assimilate(path=os.path.join(module_dir, "..", "test_files", "FF_working"), input_file="test.qin_freq1", output_file="test.qout_freq1")
        self.assertEqual(doc["input"]["job_type"],"freq")
        self.assertEqual(doc["output"]["job_type"],"freq")
        test_freqs = np.array([  12.52,   45.28,  260.96,  329.08,  531.01,  582.11,  744.91, 779.2 ,  800.47,  863.15,  928.68,  969.  , 1092.86, 1124.  ,1147.64, 1209.1 , 1387.39, 1693.97, 1913.05, 3316.2 , 3341.73])
        for ii, entry in enumerate(test_freqs):
            self.assertEqual(test_freqs[ii],doc["output"]["frequencies"][ii])
        self.assertEqual(doc["output"]["enthalpy"],37.547)
        self.assertEqual(doc["output"]["entropy"],83.81)
        self.assertEqual(doc["walltime"],394.45)
        self.assertEqual(doc["cputime"],997.39)
        self.assertEqual(doc["smiles"],"O1[C](O[Li])OC=C1")
        self.assertEqual(doc["formula_pretty"],"LiH2(CO)3")
        self.assertEqual(doc["formula_anonymous"],"AB2C3D3")
        self.assertEqual(doc["chemsys"],"C-H-Li-O")
        self.assertEqual(doc["pointgroup"],"Cs")
        self.assertIn("custodian", doc)
        self.assertIn("calcs_reversed", doc)
        self.assertIn("initial_molecule", doc["input"])
        self.assertIn("initial_molecule", doc["output"])
        self.assertIn("last_updated", doc)
        self.assertIn("dir_name", doc)
        self.assertEqual(len(doc["calcs_reversed"]),1)

    def test_assimilate_FF(self):
        drone = QChemDrone(runs=['opt_0','freq_0','opt_1','freq_1','opt_2','freq_2','opt_3','freq_3'], additional_fields={"special_run_type": "frequency_flattener"})
        doc = drone.assimilate(path=os.path.join(module_dir, "..", "test_files", "FF_working"), input_file="test.qin", output_file="test.qout")
        self.assertEqual(doc["special_run_type"],"frequency_flattener")
        self.assertEqual(doc["input"]["job_type"],"opt")
        self.assertEqual(doc["output"]["job_type"],"freq")
        test_freqs = np.array([  12.52,   45.28,  260.96,  329.08,  531.01,  582.11,  744.91, 779.2 ,  800.47,  863.15,  928.68,  969.  , 1092.86, 1124.  ,1147.64, 1209.1 , 1387.39, 1693.97, 1913.05, 3316.2 , 3341.73])
        for ii, entry in enumerate(test_freqs):
            self.assertEqual(test_freqs[ii],doc["output"]["frequencies"][ii])
        self.assertEqual(doc["output"]["enthalpy"],37.547)
        self.assertEqual(doc["output"]["entropy"],83.81)
        self.assertEqual(doc["num_frequencies_flattened"],1.0)
        self.assertEqual(doc["walltime"],935.29)
        self.assertEqual(doc["cputime"],3616.6400000000003)
        self.assertEqual(doc["smiles"],"O1[C](O[Li])OC=C1")
        self.assertEqual(doc["formula_pretty"],"LiH2(CO)3")
        self.assertEqual(doc["formula_anonymous"],"AB2C3D3")
        self.assertEqual(doc["chemsys"],"C-H-Li-O")
        self.assertEqual(doc["pointgroup"],"Cs")
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


if __name__ == "__main__":
    unittest.main()
