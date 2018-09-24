# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

from atomate.qchem.firetasks.write_inputs import WriteInputFromIOSet
from atomate.qchem.firetasks.run_calc import RunQChemCustodian
from atomate.qchem.firetasks.parse_outputs import QChemToDb
from atomate.qchem.fireworks.core import OptimizeFW, FrequencyFlatteningOptimizeFW
from atomate.utils.testing import AtomateTest
from pymatgen.io.qchem.outputs import QCOutput

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "5/23/18"
__credits__ = "Brandon Wood, Shyam Dwaraknath"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")


class TestCore(AtomateTest):
    def setUp(self, lpad=False):
        out_file = os.path.join(module_dir, "..", "..", "test_files",
                                "FF_working", "test.qout.opt_0")
        qc_out = QCOutput(filename=out_file)
        self.act_mol = qc_out.data["initial_molecule"]
        super(TestCore, self).setUp(lpad=False)

    def tearDown(self):
        pass

    def test_OptimizeFW_defaults(self):
        firework = OptimizeFW(molecule=self.act_mol)
        self.assertEqual(firework.tasks[0].as_dict(),
                         WriteInputFromIOSet(
                             molecule=self.act_mol,
                             qchem_input_set="OptSet",
                             input_file="mol.qin",
                             qchem_input_params={}).as_dict())
        self.assertEqual(firework.tasks[1].as_dict(),
                         RunQChemCustodian(
                             qchem_cmd="qchem",
                             multimode="openmp",
                             input_file="mol.qin",
                             output_file="mol.qout",
                             max_cores=32,
                             job_type="normal").as_dict())
        self.assertEqual(firework.tasks[2].as_dict(),
                         QChemToDb(
                             db_file=None,
                             input_file="mol.qin",
                             output_file="mol.qout",
                             additional_fields={
                                 "task_label": "structure optimization"
                             }).as_dict())
        self.assertEqual(firework.parents, [])
        self.assertEqual(firework.name, "structure optimization")

    def test_OptimizeFW_not_defaults(self):
        firework = OptimizeFW(
            molecule=self.act_mol,
            name="special structure optimization",
            qchem_cmd="qchem -slurm",
            multimode="mpi",
            input_file="different.qin",
            output_file="not_default.qout",
            max_cores=12,
            qchem_input_params={"pcm_dielectric": 10.0},
            db_file=os.path.join(db_dir, "db.json"),
            parents=None)
        self.assertEqual(firework.tasks[0].as_dict(),
                         WriteInputFromIOSet(
                             molecule=self.act_mol,
                             qchem_input_set="OptSet",
                             input_file="different.qin",
                             qchem_input_params={
                                 "pcm_dielectric": 10.0
                             }).as_dict())
        self.assertEqual(firework.tasks[1].as_dict(),
                         RunQChemCustodian(
                             qchem_cmd="qchem -slurm",
                             multimode="mpi",
                             input_file="different.qin",
                             output_file="not_default.qout",
                             max_cores=12,
                             job_type="normal").as_dict())
        self.assertEqual(firework.tasks[2].as_dict(),
                         QChemToDb(
                             db_file=os.path.join(db_dir, "db.json"),
                             input_file="different.qin",
                             output_file="not_default.qout",
                             additional_fields={
                                 "task_label": "special structure optimization"
                             }).as_dict())
        self.assertEqual(firework.parents, [])
        self.assertEqual(firework.name, "special structure optimization")

    def test_FrequencyFlatteningOptimizeFW_defaults(self):
        firework = FrequencyFlatteningOptimizeFW(molecule=self.act_mol)
        self.assertEqual(firework.tasks[0].as_dict(),
                         WriteInputFromIOSet(
                             molecule=self.act_mol,
                             qchem_input_set="OptSet",
                             input_file="mol.qin",
                             qchem_input_params={}).as_dict())
        self.assertEqual(firework.tasks[1].as_dict(),
                         RunQChemCustodian(
                             qchem_cmd="qchem",
                             multimode="openmp",
                             input_file="mol.qin",
                             output_file="mol.qout",
                             max_cores=32,
                             job_type="opt_with_frequency_flattener",
                             max_iterations=10,
                             max_molecule_perturb_scale=0.3,
                             reversed_direction=False).as_dict())
        self.assertEqual(firework.tasks[2].as_dict(),
                         QChemToDb(
                             db_file=None,
                             input_file="mol.qin",
                             output_file="mol.qout",
                             additional_fields={
                                 "task_label":
                                 "frequency flattening structure optimization",
                                 "special_run_type":
                                 "frequency_flattener"
                             }).as_dict())
        self.assertEqual(firework.parents, [])
        self.assertEqual(firework.name,
                         "frequency flattening structure optimization")

    def test_FrequencyFlatteningOptimizeFW_not_defaults(self):
        firework = FrequencyFlatteningOptimizeFW(
            molecule=self.act_mol,
            name="special frequency flattening structure optimization",
            qchem_cmd="qchem -slurm",
            multimode="mpi",
            input_file="different.qin",
            output_file="not_default.qout",
            max_cores=12,
            qchem_input_params={"pcm_dielectric": 10.0},
            max_iterations=5,
            max_molecule_perturb_scale=0.2,
            reversed_direction=True,
            db_file=os.path.join(db_dir, "db.json"),
            parents=None)
        self.assertEqual(firework.tasks[0].as_dict(),
                         WriteInputFromIOSet(
                             molecule=self.act_mol,
                             qchem_input_set="OptSet",
                             input_file="different.qin",
                             qchem_input_params={
                                 "pcm_dielectric": 10.0
                             }).as_dict())
        self.assertEqual(firework.tasks[1].as_dict(),
                         RunQChemCustodian(
                             qchem_cmd="qchem -slurm",
                             multimode="mpi",
                             input_file="different.qin",
                             output_file="not_default.qout",
                             max_cores=12,
                             job_type="opt_with_frequency_flattener",
                             max_iterations=5,
                             max_molecule_perturb_scale=0.2,
                             reversed_direction=True).as_dict())
        self.assertEqual(
            firework.tasks[2].as_dict(),
            QChemToDb(
                db_file=os.path.join(db_dir, "db.json"),
                input_file="different.qin",
                output_file="not_default.qout",
                additional_fields={
                    "task_label":
                    "special frequency flattening structure optimization",
                    "special_run_type":
                    "frequency_flattener"
                }).as_dict())
        self.assertEqual(firework.parents, [])
        self.assertEqual(firework.name,
                         "special frequency flattening structure optimization")


if __name__ == "__main__":
    unittest.main()
