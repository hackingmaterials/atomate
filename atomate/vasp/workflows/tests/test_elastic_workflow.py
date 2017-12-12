# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

import numpy as np

from monty.serialization import loadfn

from fireworks import FWorker, Firework, Workflow
from fireworks.core.rocket_launcher import rapidfire

from atomate.vasp.powerups import use_fake_vasp, add_modify_incar
from atomate.vasp.workflows.base.elastic import get_wf_elastic_constant
from atomate.vasp.workflows.presets.core import wf_elastic_constant, wf_elastic_constant_minimal, get_wf
from atomate.vasp.firetasks.parse_outputs import ElasticTensorToDb
from atomate.utils.testing import AtomateTest

from pymatgen.util.testing import PymatgenTest
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen import Structure

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
        self.tf_loc = os.path.join(ref_dir, 'elastic_wf')
        self.struct_si = PymatgenTest.get_structure("Si")
        self.struct_si = SpacegroupAnalyzer(self.struct_si).get_conventional_standard_structure()
        self.opt_struct = Structure.from_file(os.path.join(self.tf_loc, '1', 'inputs', 'POSCAR'))

        # Base WF
        self.base_wf = get_wf(self.struct_si, "optimize_only.yaml",
                              params=[{"db_file": ">>db_file<<"}])
        self.base_wf.append_wf(get_wf_elastic_constant(self.struct_si, db_file='>>db_file<<',
                                                       stencils=[[0.01]]*3 + [[0.03]]*3,
                                                       copy_vasp_outputs=True),
                               self.base_wf.leaf_fw_ids)
        self.base_wf_noopt = get_wf_elastic_constant(self.opt_struct, stencils=[[0.01]]*3 + [[0.03]]*3,
                                                     copy_vasp_outputs=False, sym_reduce=False,
                                                     db_file='>>db_file<<')
        ec_incar_update = {'incar_update': {'EDIFF': 1e-6, 'ENCUT': 700}}
        self.base_wf = add_modify_incar(self.base_wf, ec_incar_update)
        self.base_wf_noopt = add_modify_incar(self.base_wf_noopt, ec_incar_update)

        # Full preset WF
        self.preset_wf = wf_elastic_constant(self.struct_si)

        # Minimal WF
        self.minimal_wf = wf_elastic_constant_minimal(self.opt_struct)
        self.minimal_wf = add_modify_incar(self.minimal_wf, ec_incar_update)

        # TOEC WF (minimal)
        self.toec_wf = wf_elastic_constant_minimal(self.struct_si, order=3) 
        self.toec_wf = add_modify_incar(self.toec_wf, ec_incar_update)
        toec_data = loadfn(os.path.join(self.tf_loc, 'toec_data.json'))
        # Rather than run entire workflow, preload the spec to test the analysis
        toec_analysis = Firework([ElasticTensorToDb(structure=self.struct_si, order=3, 
                                                    db_file=">>db_file<<")], 
                                 spec={"deformation_tasks": toec_data['deformation_tasks']})
        self.toec_analysis = Workflow([toec_analysis])
        # Check 4th order to see if constructed correctly
        self.foec_wf = wf_elastic_constant_minimal(self.struct_si, order=4)

    def _simulate_vasprun(self, wf):
        reference_dir = os.path.abspath(self.tf_loc)
        if len(wf.fws) > 6:
            si_ref_dirs = {"structure optimization": os.path.join(reference_dir, "1"),
                           "elastic deformation 0": os.path.join(reference_dir, "7"),
                           "elastic deformation 1": os.path.join(reference_dir, "6"),
                           "elastic deformation 2": os.path.join(reference_dir, "5"),
                           "elastic deformation 3": os.path.join(reference_dir, "4"),
                           "elastic deformation 4": os.path.join(reference_dir, "3"),
                           "elastic deformation 5": os.path.join(reference_dir, "2")}
        elif len(wf.fws) == 3:
            si_ref_dirs = {"elastic deformation 0": os.path.join(reference_dir, "7"),
                           "elastic deformation 1": os.path.join(reference_dir, "4")}
        return use_fake_vasp(wf, si_ref_dirs, params_to_check=["ENCUT"])

    def _check_run(self, d, mode):
        if mode not in ["structure optimization", "elastic deformation 0",
                        "elastic deformation 3", "elastic analysis", "toec analysis"]:
            raise ValueError("Invalid mode!")

        if mode not in ["elastic analysis", "toec analysis"]:
            self.assertEqual(d["formula_pretty"], "Si")
            self.assertEqual(d["formula_anonymous"], "A")
            self.assertEqual(d["nelements"], 1)
            self.assertEqual(d["state"], "successful")
        
        if mode in ["structure optimization"]:
            self.assertAlmostEqual(d["calcs_reversed"][0]["output"]["structure"]["lattice"]["a"], 5.469, 2)
            self.assertAlmostEqual(d["output"]["energy_per_atom"], -5.423, 2)

        elif mode in ["elastic deformation 0"]:
            stress = np.diag([-14.716,-5.121, -5.121])
            np.testing.assert_allclose(stress,
                    d["calcs_reversed"][0]["output"]["ionic_steps"][-1]["stress"], rtol=1e-2)

        elif mode in ["elastic deformation 3"]:
            stress = d["calcs_reversed"][0]["output"]["ionic_steps"][-1]["stress"]
            self.assertAlmostEqual(stress[2][1], -45.035, places=1)

        elif mode in ["elastic analysis"]:
            c_ij = np.array(d['elastic_tensor']['ieee_format'])
            np.testing.assert_allclose([c_ij[0, 0], c_ij[0, 1], c_ij[3, 3]],
                                       [146.25, 52.26, 75.23], rtol=1e-2)
            self.assertAlmostEqual(d['derived_properties']['k_vrh'], 83.6, places=1)

        elif mode in ["toec analysis"]:
            c_ij = np.array(d['elastic_tensor']['ieee_format'][1])
            np.testing.assert_allclose([c_ij[0, 0, 0], c_ij[0, 0, 1], c_ij[0, 1, 2]],
                                       [-779.6, -432.5, -93.3], rtol=1e-1)

    def test_wf(self):
        self.base_wf = self._simulate_vasprun(self.base_wf)
        self.base_wf_noopt = self._simulate_vasprun(self.base_wf_noopt)
        self.minimal_wf = self._simulate_vasprun(self.minimal_wf)

        self.assertEqual(len(self.base_wf.fws), 8)
        self.assertEqual(len(self.base_wf_noopt.fws), 7)
        self.assertEqual(len(self.minimal_wf.fws), 3)
        self.assertEqual(len(self.toec_wf.fws), 17)
        self.assertEqual(len(self.preset_wf.fws), 26)
        self.assertEqual(len(self.foec_wf.fws), 49)

        # check vasp parameters for ionic relaxation
        defo_vis = [fw.tasks[1]['vasp_input_set'] 
                    for fw in self.base_wf.fws if "deform" in fw.name]
        assert all([vis.user_incar_settings['NSW'] == 99 for vis in defo_vis])
        assert all([vis.user_incar_settings['IBRION'] == 2 for vis in defo_vis])
        # check preset parameters
        defo_vis = [fw.tasks[2]['vasp_input_set'] 
                    for fw in self.preset_wf.fws if "deform" in fw.name]
        assert all([vis.user_incar_settings['ENCUT'] == 700 for vis in defo_vis])
        assert all([vis.user_kpoints_settings.get('grid_density') == 7000 
                    for vis in defo_vis])

        self.lp.add_wf(self.base_wf)
        self.lp.add_wf(self.base_wf_noopt)
        self.lp.add_wf(self.toec_analysis)

        rapidfire(self.lp, fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))

        # check relaxation
        d = self.get_task_collection().find_one({"task_label": "structure optimization"})
        self._check_run(d, mode="structure optimization")
        # check two of the deformation calculations
        d = self.get_task_collection().find_one({"task_label": "elastic deformation 0"})
        self._check_run(d, mode="elastic deformation 0")
        
        d = self.get_task_collection().find_one({"task_label": "elastic deformation 3"})
        self._check_run(d, mode="elastic deformation 3")

        # check the final results
        d = self.get_task_collection(coll_name="elasticity").find_one(
                {'order': 2, "optimized_structure": {"$exists":True}})
        self._check_run(d, mode="elastic analysis")

        # check third-order results
        d = self.get_task_collection(coll_name="elasticity").find_one({'order': 3})
        self._check_run(d, mode="toec analysis")

        wf = self.lp.get_wf_by_fw_id(1)
        self.assertTrue(all([s == 'COMPLETED' for s in wf.fw_states.values()]))


if __name__ == "__main__":
    unittest.main()
