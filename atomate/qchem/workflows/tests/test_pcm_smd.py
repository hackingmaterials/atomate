# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

from fireworks import FWorker
from fireworks.core.rocket_launcher import rapidfire
from atomate.utils.testing import AtomateTest
from pymatgen.core import Molecule
from pymatgen.io.qchem.inputs import QCInput
from atomate.qchem.powerups import use_fake_qchem
from atomate.qchem.workflows.base.pcm_smd import get_wf_pcm_smd

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "4/22/19"


module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")


class TestPCMSMD(AtomateTest):
    def test_pcm_smd(self):
        # location of test files
        test_double_FF_files = os.path.join(module_dir, "..", "..",
                                            "test_files", "double_FF_wf")
        # define starting molecule and workflow object
        initial_qcin = QCInput.from_file(
            os.path.join(test_double_FF_files, "block", "launcher_first",
                         "mol.qin.opt_0"))
        initial_mol = initial_qcin.molecule

        try:
            wf = get_wf_pcm_smd(
                molecule=initial_mol,
                pcm_dielectric=10.0,
                smd_solvent="custom",
                qchem_input_params={
                    "basis_set": "6-311++g**",
                    "scf_algorithm": "diis",
                    "overwrite_inputs": {
                        "rem": {
                            "sym_ignore": "true"
                        }
                    }
                })
            self.assertEqual(True,False)
        except RuntimeError:
            wf = get_wf_pcm_smd(
                molecule=initial_mol,
                pcm_dielectric=90.0,
                smd_solvent="custom",
                qchem_input_params={
                    "basis_set": "6-311++g**",
                    "scf_algorithm": "diis",
                    "custom_smd": "90.00,1.415,0.00,0.735,20.2,0.00,0.00",
                    "overwrite_inputs": {
                        "rem": {
                            "sym_ignore": "true"
                        }
                    }
                })

if __name__ == "__main__":
    unittest.main()
