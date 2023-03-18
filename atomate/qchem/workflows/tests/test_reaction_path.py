import copy
import os
import unittest

import numpy as np
from fireworks import FWorker
from fireworks.core.rocket_launcher import rapidfire
from pymatgen.core import Molecule
from pymatgen.io.qchem.outputs import QCOutput

from atomate.qchem.powerups import use_fake_qchem
from atomate.qchem.workflows.base.reaction_path import get_wf_reaction_path_with_ts
from atomate.utils.testing import AtomateTest

__author__ = "Evan Spotte-Smith"
__copyright__ = "Copyright 2021, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Alpha"
__date__ = "02/04/2021"


module_dir = os.path.dirname(os.path.abspath(__file__))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")


class TestReactionPath(AtomateTest):
    def test_reaction_path_with_ts(self):
        # location of test files
        test_double_FF_files = os.path.join(
            module_dir, "..", "..", "test_files", "reaction_path_wf"
        )
        # define starting molecule and workflow object
        initial_qout = QCOutput(os.path.join(test_double_FF_files, "..", "ts.out"))
        initial_mol = initial_qout.data["initial_molecule"]
        formula = initial_mol.composition.alphabetical_formula
        mode = initial_qout.data["frequency_mode_vectors"][0]

        real_wf = get_wf_reaction_path_with_ts(
            molecule=initial_mol,
            mode=mode,
            scale=0.6,
            suffix="HERE",
            qchem_input_params={
                "dft_rung": 4,
                "basis_set": "def2-tzvppd",
                "smd_solvent": "custom",
                "custom_smd": "18.5,1.415,0.00,0.735,20.2,0.00,0.00",
                "overwrite_inputs": {
                    "rem": {
                        "scf_algorithm": "diis",
                        "thresh": 14,
                        "method": "wb97xv",
                    }
                },
            },
        )
        # use powerup to replace run with fake run
        ref_dirs = {
            f"{formula}:perturb_forwardsHERE": os.path.join(
                test_double_FF_files, "block", "launcher_forwards"
            ),
            f"{formula}:perturb_backwardsHERE": os.path.join(
                test_double_FF_files, "block", "launcher_backwards"
            ),
        }
        fake_wf = use_fake_qchem(real_wf, ref_dirs)
        self.lp.add_wf(fake_wf)
        rapidfire(
            self.lp,
            fworker=FWorker(
                env={"max_cores": 32, "db_file": os.path.join(db_dir, "db.json")}
            ),
        )

        wf_test = self.lp.get_wf_by_fw_id(1)
        self.assertTrue(all([s == "COMPLETED" for s in wf_test.fw_states.values()]))

        forwards = self.get_task_collection().find_one(
            {"task_label": f"{formula}:perturb_forwardsHERE"}
        )
        forwards_final_mol = Molecule.from_dict(forwards["input"]["initial_molecule"])

        backwards = self.get_task_collection().find_one(
            {"task_label": f"{formula}:perturb_backwardsHERE"}
        )
        backwards_final_mol = Molecule.from_dict(backwards["input"]["initial_molecule"])

        mol_copy_for = copy.deepcopy(initial_mol)
        mol_copy_back = copy.deepcopy(initial_mol)

        for ii in range(len(mol_copy_for)):
            vec = np.array(mode[ii])
            mol_copy_for.translate_sites(indices=[ii], vector=vec * 0.6)

        for ii in range(len(mol_copy_back)):
            vec = np.array(mode[ii])
            mol_copy_back.translate_sites(indices=[ii], vector=vec * -0.6)

        np.testing.assert_allclose(
            mol_copy_for.cart_coords, forwards_final_mol.cart_coords, atol=0.0001
        )

        np.testing.assert_allclose(
            mol_copy_back.cart_coords, backwards_final_mol.cart_coords, atol=0.0001
        )


if __name__ == "__main__":
    unittest.main()
