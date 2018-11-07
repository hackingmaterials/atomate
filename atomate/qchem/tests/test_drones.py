# coding: utf-8
# Copyright (c) Materials Virtual Lab.
# Distributed under the terms of the BSD License.

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
from atomate.qchem.drones import QChemDrone
from pymatgen.core.structure import Molecule
import numpy as np
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.analysis.graphs import MoleculeGraph

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "4/29/18"
__credits__ = "Brandon Wood, Shyam Dwaraknath"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class QChemDroneTest(unittest.TestCase):
    def test_assimilate_opt(self):
        drone = QChemDrone()
        doc = drone.assimilate(
            path=os.path.join(module_dir, "..", "test_files", "FF_working"),
            input_file="test.qin.opt_1",
            output_file="test.qout.opt_1",
            multirun=False)
        self.assertEqual(doc["input"]["job_type"], "opt")
        self.assertEqual(doc["output"]["job_type"], "opt")
        self.assertEqual(doc["output"]["final_energy"], -348.652462579636)
        self.assertEqual(doc["walltime"], 62.83)
        self.assertEqual(doc["cputime"], 715.76)
        self.assertEqual(doc["smiles"], "O1[C](O[Li])OC=C1")
        self.assertEqual(doc["formula_pretty"], "LiH2(CO)3")
        self.assertEqual(doc["formula_anonymous"], "AB2C3D3")
        self.assertEqual(doc["chemsys"], "C-H-Li-O")
        self.assertEqual(doc["pointgroup"], "Cs")
        self.assertIn("custodian", doc)
        self.assertIn("calcs_reversed", doc)
        self.assertIn("initial_molecule", doc["input"])
        self.assertIn("initial_molecule", doc["output"])
        self.assertIn("optimized_molecule", doc["output"])
        self.assertIn("last_updated", doc)
        self.assertIn("dir_name", doc)
        self.assertEqual(len(doc["calcs_reversed"]), 1)

    def test_assimilate_freq(self):
        drone = QChemDrone()
        doc = drone.assimilate(
            path=os.path.join(module_dir, "..", "test_files", "FF_working"),
            input_file="test.qin.freq_1",
            output_file="test.qout.freq_1",
            multirun=False)
        self.assertEqual(doc["input"]["job_type"], "freq")
        self.assertEqual(doc["output"]["job_type"], "freq")
        test_freqs = np.array([
            12.52, 45.28, 260.96, 329.08, 531.01, 582.11, 744.91, 779.2,
            800.47, 863.15, 928.68, 969., 1092.86, 1124., 1147.64, 1209.1,
            1387.39, 1693.97, 1913.05, 3316.2, 3341.73
        ])
        for ii in enumerate(test_freqs):
            self.assertEqual(test_freqs[ii[0]], doc["output"]["frequencies"][ii[0]])
        self.assertEqual(doc["output"]["enthalpy"], 37.547)
        self.assertEqual(doc["output"]["entropy"], 83.81)
        self.assertEqual(doc["walltime"], 394.45)
        self.assertEqual(doc["cputime"], 997.39)
        self.assertEqual(doc["smiles"], "O1[C](O[Li])OC=C1")
        self.assertEqual(doc["formula_pretty"], "LiH2(CO)3")
        self.assertEqual(doc["formula_anonymous"], "AB2C3D3")
        self.assertEqual(doc["chemsys"], "C-H-Li-O")
        self.assertEqual(doc["pointgroup"], "Cs")
        self.assertIn("custodian", doc)
        self.assertIn("calcs_reversed", doc)
        self.assertIn("initial_molecule", doc["input"])
        self.assertIn("initial_molecule", doc["output"])
        self.assertIn("last_updated", doc)
        self.assertIn("dir_name", doc)
        self.assertEqual(len(doc["calcs_reversed"]), 1)

    def test_assimilate_FF(self):
        drone = QChemDrone(
            runs=[
                "opt_0", "freq_0", "opt_1", "freq_1", "opt_2", "freq_2",
                "opt_3", "freq_3"
            ],
            additional_fields={"special_run_type": "frequency_flattener"})
        doc = drone.assimilate(
            path=os.path.join(module_dir, "..", "test_files", "FF_working"),
            input_file="test.qin",
            output_file="test.qout",
            multirun=False)
        self.assertEqual(doc["special_run_type"], "frequency_flattener")
        self.assertEqual(doc["input"]["job_type"], "opt")
        self.assertEqual(doc["output"]["job_type"], "freq")
        test_freqs = np.array([
            12.52, 45.28, 260.96, 329.08, 531.01, 582.11, 744.91, 779.2,
            800.47, 863.15, 928.68, 969., 1092.86, 1124., 1147.64, 1209.1,
            1387.39, 1693.97, 1913.05, 3316.2, 3341.73
        ])
        for ii in enumerate(test_freqs):
            self.assertEqual(test_freqs[ii[0]], doc["output"]["frequencies"][ii[0]])
            self.assertEqual(doc["output"]["frequencies"][ii[0]],
                             doc["calcs_reversed"][0]["frequencies"][ii[0]])
        self.assertEqual(doc["output"]["enthalpy"], 37.547)
        self.assertEqual(doc["output"]["entropy"], 83.81)
        self.assertEqual(doc["num_frequencies_flattened"], 1)
        self.assertEqual(doc["walltime"], 935.29)
        self.assertEqual(doc["cputime"], 3616.6400000000003)
        self.assertEqual(doc["smiles"], "O1[C](O[Li])OC=C1")
        self.assertEqual(doc["formula_pretty"], "LiH2(CO)3")
        self.assertEqual(doc["formula_anonymous"], "AB2C3D3")
        self.assertEqual(doc["chemsys"], "C-H-Li-O")
        self.assertEqual(doc["pointgroup"], "Cs")
        self.assertIn("custodian", doc)
        self.assertIn("calcs_reversed", doc)
        self.assertIn("initial_molecule", doc["input"])
        self.assertIn("initial_molecule", doc["output"])
        self.assertIn("optimized_molecule", doc["output"])
        self.assertIn("last_updated", doc)
        self.assertIn("dir_name", doc)
        self.assertEqual(len(doc["calcs_reversed"]), 4)
        self.assertEqual(
            list(doc["calcs_reversed"][0].keys()),
            list(doc["calcs_reversed"][2].keys()))
        self.assertEqual(
            list(doc["calcs_reversed"][1].keys()),
            list(doc["calcs_reversed"][3].keys()))

    def test_assimilate_bad_FF(self):
        drone = QChemDrone(additional_fields={"special_run_type": "frequency_flattener"})
        doc = drone.assimilate(
            path=os.path.join(module_dir, "..", "test_files", "launcher_bad_FF"),
            input_file="mol.qin",
            output_file="mol.qout",
            multirun=False)
        self.assertEqual(doc["special_run_type"], "frequency_flattener")
        self.assertEqual(doc["input"]["job_type"], "opt")
        self.assertEqual(doc["output"]["job_type"], "freq")
        self.assertEqual(doc["state"], "unsuccessful")

    def test_multirun(self):
        drone = QChemDrone()
        doc = drone.assimilate(
            path=os.path.join(module_dir, "..", "test_files", "julian_nt"),
            input_file="julian.qin",
            output_file="julian.qout",
            multirun=True)
        self.assertEqual(doc["input"]["job_type"], "optimization")
        self.assertEqual(doc["output"]["job_type"], "frequency")
        test_freqs = np.array([
            -69.17, 117.81, 244.67, 257.93, 530., 579.64, 737.42, 771.1,
            787.32, 869.29, 924.77, 962.67, 1084.55, 1117.49, 1143.1, 1196.27,
            1378.76, 1696.26, 1860.75, 3321.43
        ])
        for ii in enumerate(test_freqs):
            self.assertEqual(test_freqs[ii[0]], doc["output"]["frequencies"][ii[0]])
            self.assertEqual(doc["output"]["frequencies"][ii[0]],
                             doc["calcs_reversed"][0]["frequencies"][ii[0]])
        self.assertEqual(doc["output"]["enthalpy"], 36.755)
        self.assertEqual(doc["output"]["entropy"], 74.989)
        self.assertEqual(doc["walltime"], 684.6300000000001)
        self.assertEqual(doc["cputime"], 4039.37)
        self.assertEqual(doc["smiles"], "O1[C](O[Li])OC=C1")
        self.assertEqual(doc["formula_pretty"], "LiH2(CO)3")
        self.assertEqual(doc["formula_anonymous"], "AB2C3D3")
        self.assertEqual(doc["chemsys"], "C-H-Li-O")
        self.assertEqual(doc["pointgroup"], "C2")
        self.assertIn("calcs_reversed", doc)
        self.assertIn("initial_molecule", doc["input"])
        self.assertIn("initial_molecule", doc["output"])
        self.assertIn("optimized_molecule", doc["output"])
        self.assertIn("last_updated", doc)
        self.assertIn("dir_name", doc)
        self.assertEqual(len(doc["calcs_reversed"]), 3)
        self.assertEqual(doc["calcs_reversed"][0]["task"]["name"], "calc2")
        self.assertEqual(doc["calcs_reversed"][-1]["task"]["name"], "calc0")

    def test_assimilate_unstable_opt(self):
        drone = QChemDrone(
            runs=[
                "opt_0", "freq_0", "opt_1", "freq_1", "opt_2", "freq_2",
                "opt_3", "freq_3"
            ],
            additional_fields={"special_run_type": "frequency_flattener"})
        doc = drone.assimilate(
            path=os.path.join(module_dir, "..", "test_files", "2620_complete"),
            input_file="mol.qin",
            output_file="mol.qout",
            multirun=False)
        self.assertEqual(doc["input"]["job_type"], "opt")
        self.assertEqual(doc["output"]["job_type"], "opt")
        self.assertEqual(doc["output"]["final_energy"], "unstable")
        self.assertEqual(doc["smiles"], "[S](=O)[N]S[C]")
        self.assertEqual(doc["state"], "unsuccessful")
        self.assertEqual(doc["num_frequencies_flattened"], 0)
        self.assertEqual(doc["walltime"], None)
        self.assertEqual(doc["cputime"], None)
        self.assertEqual(doc["formula_pretty"], "CS2NO")
        self.assertEqual(doc["formula_anonymous"], "ABCD2")
        self.assertEqual(doc["chemsys"], "C-N-O-S")
        self.assertEqual(doc["pointgroup"], "C1")
        self.assertEqual(doc["orig"]["rem"], doc["calcs_reversed"][-1]["input"]["rem"])
        self.assertEqual(doc["orig"]["molecule"], doc["calcs_reversed"][-1]["input"]["molecule"])
        orig_molgraph = MoleculeGraph.with_local_env_strategy(Molecule.from_dict(doc["orig"]["molecule"]),
                                                              OpenBabelNN(),
                                                              reorder=False,
                                                              extend_structure=False)
        initial_molgraph = MoleculeGraph.with_local_env_strategy(Molecule.from_dict(doc["input"]["initial_molecule"]),
                                                                 OpenBabelNN(),
                                                                 reorder=False,
                                                                 extend_structure=False)
        self.assertEqual(orig_molgraph.isomorphic_to(initial_molgraph), True)

    def test_assimilate_opt_with_hidden_changes_from_handler(self):
        drone = QChemDrone(additional_fields={"special_run_type": "frequency_flattener"})
        doc = drone.assimilate(
            path=os.path.join(module_dir, "..", "test_files", "1746_complete"),
            input_file="mol.qin",
            output_file="mol.qout",
            multirun=False)
        self.assertEqual(doc["input"]["job_type"], "opt")
        self.assertEqual(doc["output"]["job_type"], "freq")
        self.assertEqual(doc["output"]["final_energy"], -303.835532370106)
        self.assertEqual(doc["smiles"], "O1C(=CC1=O)[CH]")
        self.assertEqual(doc["state"], "successful")
        self.assertEqual(doc["num_frequencies_flattened"], 0)
        self.assertEqual(doc["walltime"], 631.54)
        self.assertEqual(doc["cputime"], 7471.17)
        self.assertEqual(doc["formula_pretty"], "HC2O")
        self.assertEqual(doc["formula_anonymous"], "ABC2")
        self.assertEqual(doc["chemsys"], "C-H-O")
        self.assertEqual(doc["pointgroup"], "C1")
        self.assertEqual(doc["orig"]["rem"], doc["calcs_reversed"][-1]["input"]["rem"])
        orig_molgraph = MoleculeGraph.with_local_env_strategy(Molecule.from_dict(doc["orig"]["molecule"]),
                                                              OpenBabelNN(),
                                                              reorder=False,
                                                              extend_structure=False)
        initial_molgraph = MoleculeGraph.with_local_env_strategy(Molecule.from_dict(doc["input"]["initial_molecule"]),
                                                                 OpenBabelNN(),
                                                                 reorder=False,
                                                                 extend_structure=False)
        self.assertEqual(orig_molgraph.isomorphic_to(initial_molgraph), False)

    def test_assimilate_disconnected_opt(self):
        drone = QChemDrone(additional_fields={"special_run_type": "frequency_flattener"})
        doc = drone.assimilate(
            path=os.path.join(module_dir, "..", "test_files", "disconnected_but_converged"),
            input_file="mol.qin",
            output_file="mol.qout",
            multirun=False)
        self.assertEqual(doc["input"]["job_type"], "opt")
        self.assertEqual(doc["output"]["job_type"], "freq")
        self.assertEqual(doc["output"]["final_energy"], -303.07602688705)
        self.assertEqual(doc["smiles"], "O=C.O=C=O")
        self.assertEqual(doc["state"], "successful")
        self.assertEqual(doc["num_frequencies_flattened"], 0)
        self.assertEqual(doc["walltime"], 492.42999999999995)
        self.assertEqual(doc["cputime"], 8825.76)
        self.assertEqual(doc["formula_pretty"], "H2C2O3")
        self.assertEqual(doc["formula_anonymous"], "A2B2C3")
        self.assertEqual(doc["chemsys"], "C-H-O")
        self.assertEqual(doc["pointgroup"], "C1")
        self.assertEqual(doc["orig"]["rem"], doc["calcs_reversed"][-1]["input"]["rem"])
        self.assertEqual(doc["calcs_reversed"][-1]["structure_change"], "unconnected_fragments")

    def test_assimilate_sp(self):
        drone = QChemDrone()
        doc = drone.assimilate(
            path=os.path.join(module_dir, "..", "test_files", "launcher_sp"),
            input_file="mol.qin",
            output_file="mol.qout",
            multirun=False)
        self.assertEqual(doc["input"]["job_type"], "sp")
        self.assertEqual(doc["output"]["job_type"], "sp")
        self.assertEqual(doc["output"]["final_energy"], -75.1151765884)
        self.assertEqual(doc["walltime"], 4.69)
        self.assertEqual(doc["cputime"], 134.03)
        self.assertEqual(doc["smiles"], "[O]")
        self.assertEqual(doc["formula_pretty"], "O2")
        self.assertEqual(doc["formula_anonymous"], "A")
        self.assertEqual(doc["chemsys"], "O")
        self.assertEqual(doc["pointgroup"], "Kh")
        self.assertIn("custodian", doc)
        self.assertIn("calcs_reversed", doc)
        self.assertIn("initial_molecule", doc["input"])
        self.assertIn("initial_molecule", doc["output"])
        self.assertIn("last_updated", doc)
        self.assertIn("dir_name", doc)
        self.assertEqual(len(doc["calcs_reversed"]), 1)

    def test_sp_with_orig(self):
        drone = QChemDrone()
        doc = drone.assimilate(
            path=os.path.join(module_dir, "..", "test_files", "launcher_bad_sp"),
            input_file="mol.qin",
            output_file="mol.qout",
            multirun=False)
        self.assertEqual(doc["input"]["job_type"], "sp")
        self.assertEqual(doc["output"]["job_type"], "sp")
        self.assertEqual(doc["output"]["final_energy"], -74.540726551)
        self.assertEqual(doc["walltime"], 3.9)
        self.assertEqual(doc["cputime"], 113.27)
        self.assertEqual(doc["smiles"], "[O]")
        self.assertEqual(doc["formula_pretty"], "O2")
        self.assertEqual(doc["formula_anonymous"], "A")
        self.assertEqual(doc["chemsys"], "O")
        self.assertEqual(doc["pointgroup"], "Kh")
        self.assertIn("custodian", doc)
        self.assertIn("calcs_reversed", doc)
        self.assertIn("initial_molecule", doc["input"])
        self.assertIn("initial_molecule", doc["output"])
        self.assertIn("last_updated", doc)
        self.assertIn("dir_name", doc)
        self.assertEqual(len(doc["calcs_reversed"]), 1)


if __name__ == "__main__":
    unittest.main()
