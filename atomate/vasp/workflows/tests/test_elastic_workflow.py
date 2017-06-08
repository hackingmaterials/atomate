# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

import numpy as np

from fireworks import FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.vasp.powerups import use_fake_vasp, add_modify_incar
from atomate.vasp.workflows.presets.core import wf_elastic_constant
from atomate.vasp.workflows.base.elastic import get_wf_elastic_constant
from atomate.utils.testing import AtomateTest

from pymatgen.util.testing import PymatgenTest
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

__author__ = 'Kiran Mathew, Joseph Montoya'
__email__ = 'montoyjh@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
ref_dir = os.path.join(module_dir, "..", "..", "test_files")

DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...


class TestElasticWorkflow(AtomateTest):

    def setUp(self):
        super(TestElasticWorkflow, self).setUp()
        self.struct_si = SpacegroupAnalyzer(PymatgenTest.get_structure("Si")).get_conventional_standard_structure()
        self.elastic_config = {"NORM_DEFORMATIONS":[0.01],
                               "SHEAR_DEFORMATIONS":[0.03],
                               "VASP_CMD": ">>vasp_cmd<<",
                               "DB_FILE": ">>db_file<<"}
        self.wf = wf_elastic_constant(self.struct_si, self.elastic_config)
        self.wf_noopt = get_wf_elastic_constant(self.struct_si, norm_deformations=[0.01],
                                                shear_deformations=[0.03], copy_vasp_outputs=False)
        mip = {"incar_update": {"ENCUT": 700}}
        self.wf_noopt = add_modify_incar(self.wf_noopt, modify_incar_params=mip)

    def _simulate_vasprun(self, wf):
        reference_dir = os.path.abspath(os.path.join(ref_dir, "elastic_wf"))
        si_ref_dirs = {"structure optimization": os.path.join(reference_dir, "1"),
                       "elastic deformation 0": os.path.join(reference_dir, "7"),
                       "elastic deformation 1": os.path.join(reference_dir, "6"),
                       "elastic deformation 2": os.path.join(reference_dir, "5"),
                       "elastic deformation 3": os.path.join(reference_dir, "4"),
                       "elastic deformation 4": os.path.join(reference_dir, "3"),
                       "elastic deformation 5": os.path.join(reference_dir, "2")}
        return use_fake_vasp(wf, si_ref_dirs, params_to_check=["ENCUT"])

    def _check_run(self, d, mode):
        if mode not in ["structure optimization", "elastic deformation 0",
                        "elastic deformation 3", "elastic analysis"]:
            raise ValueError("Invalid mode!")

        if mode not in ["elastic analysis"]:
            self.assertEqual(d["formula_pretty"], "Si")
            self.assertEqual(d["formula_anonymous"], "A")
            self.assertEqual(d["nelements"], 1)
            self.assertEqual(d["state"], "successful")
        
        if mode in ["structure optimization"]:
            self.assertAlmostEqual(d["calcs_reversed"][0]["output"]["structure"]["lattice"]["a"], 5.469, 2)
            self.assertAlmostEqual(d["output"]["energy_per_atom"], -5.425, 2)

        elif mode in ["elastic deformation 0"]:
            stress = np.diag([-14.741,-5.107, -5.107])
            np.testing.assert_allclose(stress,
                    d["calcs_reversed"][0]["output"]["ionic_steps"][-1]["stress"], rtol=1e-2)

        elif mode in ["elastic deformation 3"]:
            stress = d["calcs_reversed"][0]["output"]["ionic_steps"][-1]["stress"]
            self.assertAlmostEqual(stress[0][1], -22.4, places=1)

        elif mode in ["elastic analysis"]:
            c_ij = np.array(d['elastic_tensor'])
            np.testing.assert_allclose([c_ij[0, 0], c_ij[0, 1], c_ij[3, 3]],
                                       [146.68, 50.817, 74.706], rtol=1e-2)
            self.assertAlmostEqual(d['k_voigt'], 83, places=0)

    def test_wf(self):
        self.wf = self._simulate_vasprun(self.wf)
        self.wf_noopt = self._simulate_vasprun(self.wf_noopt)

        self.assertEqual(len(self.wf.fws), 8)
        # check vasp parameters for ionic relaxation
        defo_vis = [fw.tasks[2]['vasp_input_set']
                    for fw in self.wf.fws if "deform" in fw.name]
        assert all([vis.user_incar_settings['NSW'] == 99 for vis in defo_vis])
        assert all([vis.user_incar_settings['IBRION'] == 2 for vis in defo_vis])
        self.lp.add_wf(self.wf)
        self.lp.add_wf(self.wf_noopt)
        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))

        # check relaxation
        d = self.get_task_collection().find_one({"task_label": "elastic structure optimization"})
        self._check_run(d, mode="structure optimization")
        # check two of the deformation calculations
        d = self.get_task_collection().find_one({"task_label": "elastic deformation 0"})
        self._check_run(d, mode="elastic deformation 0")
        
        d = self.get_task_collection().find_one({"task_label": "elastic deformation 3"})
        self._check_run(d, mode="elastic deformation 3")

        # check the final results
        d = self.get_task_collection(coll_name="elasticity").find_one()
        self._check_run(d, mode="elastic analysis")

        wf = self.lp.get_wf_by_fw_id(1)
        self.assertTrue(all([s == 'COMPLETED' for s in wf.fw_states.values()]))


if __name__ == "__main__":
    unittest.main()
