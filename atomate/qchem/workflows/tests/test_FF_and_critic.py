# coding: utf-8


import os
import unittest
from monty.serialization import loadfn
from fireworks import FWorker
from fireworks.core.rocket_launcher import rapidfire
from atomate.utils.testing import AtomateTest
from pymatgen.core import Molecule
from pymatgen.io.qchem.inputs import QCInput
from atomate.qchem.powerups import use_fake_qchem
from atomate.qchem.workflows.base.FF_and_critic import get_wf_FFopt_and_critic

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "11/20/19"


module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")


class TestFFOptandCritic(AtomateTest):
    def test_FFopt_and_critic(self):
        # location of test files
        test_files = os.path.join(module_dir, "..", "..", "test_files", "critic_test_files")
        # define starting molecule and workflow object
        initial_qcin = QCInput.from_file(
            os.path.join(test_files, "FFopt", "mol.qin.orig"))
        initial_mol = initial_qcin.molecule

        real_wf = get_wf_FFopt_and_critic(
            molecule=initial_mol,
            suffix="testing",
            qchem_input_params={
                "dft_rung": 4,
                "smd_solvent": "custom",
                "custom_smd": "18.5,1.415,0.00,0.735,20.2,0.00,0.00",
                "overwrite_inputs": {
                    "rem": {
                        "thresh": "14",
                        "scf_guess_always": "True"
                    }
                }
            })
        # use powerup to replace run with fake run
        ref_dirs = {
            "{}:{}".format(initial_mol.composition.alphabetical_formula, "FFopt_testing"):
            os.path.join(test_files, "FFopt"),
            "{}:{}".format(initial_mol.composition.alphabetical_formula, "CC2_testing"):
            os.path.join(test_files, "critic_example")
        }
        fake_wf = use_fake_qchem(real_wf, ref_dirs)
        self.lp.add_wf(fake_wf)
        rapidfire(
            self.lp,
            fworker=FWorker(env={"max_cores": 32, "db_file": os.path.join(db_dir, "db.json")}))

        wf_test = self.lp.get_wf_by_fw_id(1)
        self.assertTrue(
            all([s == "COMPLETED" for s in wf_test.fw_states.values()]))

        FFopt = self.get_task_collection().find_one({
            "task_label":
            "{}:{}".format(initial_mol.composition.alphabetical_formula, "FFopt_testing")
        })
        self.assertEqual(FFopt["calcs_reversed"][0]["input"]["smx"]["solvent"],
                         "other")
        self.assertEqual(FFopt["num_frequencies_flattened"], 0)
        FFopt_final_mol = Molecule.from_dict(
            FFopt["output"]["optimized_molecule"])

        CC2 = self.get_task_collection().find_one({
            "task_label":
            "{}:{}".format(initial_mol.composition.alphabetical_formula, "CC2_testing")
        })
        CC2_initial_mol = Molecule.from_dict(
            CC2["input"]["initial_molecule"])

        self.assertEqual(FFopt_final_mol, CC2_initial_mol)
        self.assertEqual(CC2["output"]["job_type"], "sp")
        self.assertEqual(CC2["output"]["final_energy"], -343.4820411597)
        critic2_drone_ref = loadfn(os.path.join(test_files, "critic_example", "critic2_drone_ref.json"))
        self.assertEqual(CC2["critic2"], critic2_drone_ref)


if __name__ == "__main__":
    unittest.main()
