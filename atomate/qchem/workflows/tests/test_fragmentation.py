# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
import json

from fireworks import FWorker
from fireworks.core.rocket_launcher import rapidfire
from atomate.utils.testing import AtomateTest
from pymatgen.io.qchem.inputs import QCInput
from atomate.qchem.powerups import use_fake_qchem
from atomate.qchem.workflows.base.fragmentation import get_fragmentation_wf
from atomate.qchem.database import QChemCalcDb
from pymatgen.core import Molecule
try:
    from unittest.mock import patch, MagicMock
except ImportError:
    from mock import patch, MagicMock

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "6/1/18"
__credits__ = "Brandon Wood, Shyam Dwaraknath"


module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")


class TestFragmentation(AtomateTest):
    def test_Fragmentation(self):
        with patch("atomate.qchem.firetasks.fragmenter.FWAction") as FWAction_patch:
            mock_FWAction = MagicMock()
            FWAction_patch.return_value = mock_FWAction
            mock_FWAction.as_dict.return_value = {'stored_data': {}, 'exit': False, 'update_spec': {}, 'mod_spec': [], 'additions': [], 'detours': [], 'defuse_children': False, 'defuse_workflow': False}

            # location of test files
            test_FF_then_fragment_files = os.path.join(module_dir, "..", "..",
                                                "test_files", "FF_then_fragment_wf")
            # define starting molecule and workflow object
            initial_qcin = QCInput.from_file(
                os.path.join(test_FF_then_fragment_files, "block", "launcher_first",
                             "mol.qin.opt_0"))
            initial_mol = initial_qcin.molecule
            real_wf = get_fragmentation_wf(molecule=initial_mol, depth=0, do_triplets=False)
            # use powerup to replace run with fake run
            ref_dirs = {
                "first FF":
                os.path.join(test_FF_then_fragment_files, "block", "launcher_first"),
                "fragment and FF_opt":
                os.path.join(test_FF_then_fragment_files, "block", "launcher_second")
            }
            fake_wf = use_fake_qchem(real_wf, ref_dirs)
            self.lp.add_wf(fake_wf)
            rapidfire(
                self.lp,
                fworker=FWorker(env={"max_cores": 32, "db_file": os.path.join(db_dir, "db.json")}), pdb_on_exception=True)

            first_FF = self.get_task_collection().find_one({
                "task_label":
                "first FF"
            })
            self.assertEqual(first_FF["calcs_reversed"][0]["input"]["solvent"],
                             None)
            self.assertEqual(first_FF["num_frequencies_flattened"], 0)
            self.assertEqual(len(FWAction_patch.call_args[1]["additions"]), 5 * 3)

    def test_no_opt_Fragmentation(self):
        db_file = os.path.join(db_dir, "db.json")
        mmdb = QChemCalcDb.from_db_file(db_file, admin=True)
        with open(os.path.join(module_dir, "..", "..", "test_files","sb40.json")) as f:
            tmp = json.load(f)
            for entry in tmp:
                mmdb.insert(entry)
        with patch("atomate.qchem.firetasks.fragmenter.FWAction") as FWAction_patch:
            mock_FWAction = MagicMock()
            FWAction_patch.return_value = mock_FWAction
            mock_FWAction.as_dict.return_value = {'stored_data': {}, 'exit': False, 'update_spec': {}, 'mod_spec': [], 'additions': [], 'detours': [], 'defuse_children': False, 'defuse_workflow': False}

            # define starting molecule and workflow object
            initial_mol = Molecule.from_file(os.path.join(module_dir, "..", "..", "test_files", "top_11", "EC.xyz"))
            initial_mol.set_charge_and_spin(charge=-1)
            wf = get_fragmentation_wf(molecule=initial_mol, depth=1, pcm_dielectric=40.0, do_optimization=False, check_db=True)
            self.lp.add_wf(wf)
            rapidfire(
                self.lp,
                fworker=FWorker(env={"max_cores": 24, "db_file": db_file}), pdb_on_exception=True)

            self.assertEqual(len(FWAction_patch.call_args[1]["additions"]), 0)
        mmdb.reset()


if __name__ == "__main__":
    unittest.main()
