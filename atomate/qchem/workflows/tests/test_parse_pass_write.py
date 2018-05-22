# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

from atomate.qchem.firetasks.geo_transformations import RotateTorsion
from atomate.qchem.firetasks.write_inputs import WriteInputFromIOSet
from atomate.qchem.firetasks.parse_outputs import QChemToDb
from fireworks import Firework, Workflow, FWorker
from fireworks.core.rocket_launcher import rapidfire
from atomate.utils.testing import AtomateTest
from pymatgen.core import Molecule
import numpy as np

__author__ = 'Brandon Wood'
__email__ = 'b.wood@berkeley.edu'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")


class TestParsePassWrite(AtomateTest):
    @classmethod
    def setUpClass(cls):

        cls.act_mol = Molecule.from_file(os.path.join(module_dir, "..", "..", "test_files",
                                                      "FF_working", "test.qin.opt_1"))

    def tearDown(self):
        # this removes the scratch dir made by AtomateTest
        shutil.rmtree(self.scratch_dir)
        # this removes the file that gets written
        for x in ["mol.qin"]:
            if os.path.exists(os.path.join(module_dir, x)):
                os.remove(os.path.join(module_dir, x))

    def test_parse_pass_write(self):

        input = "test.qin.opt_1"
        output = "test.qout.opt_1"
        calc_dir = os.path.join(module_dir, "..", "..", "test_files", "FF_working")

        p_task = QChemToDb(calc_dir=calc_dir, input_file=input, output_file=output)
        fw1 = Firework([p_task])
        w_task = WriteInputFromIOSet(qchem_input_set="OptSet", write_to_dir=module_dir)
        fw2 = Firework([w_task], parents=fw1)
        wf = Workflow([fw1, fw2])

        self.lp.add_wf(wf)
        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))

        test_mol = QCInput.from_file(os.path.join(module_dir, "mol.qin"))
        np.testing.assert_equal(self.act_mol.species, test_mol.species)
        np.testing.assert_equal(self.act_mol.cart_coords, test_mol.cart_coords)

if __name__ == '__main__':
    unittest.main()
