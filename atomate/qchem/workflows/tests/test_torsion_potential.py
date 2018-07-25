# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

from fireworks import FWorker
from fireworks.core.rocket_launcher import rapidfire
from atomate.utils.testing import AtomateTest
from pymatgen.core import Molecule
from pymatgen.io.qchem_io.outputs import QCOutput
from atomate.qchem.workflows.base.torsion_potential import get_wf_torsion_potential
from atomate.qchem.powerups import use_fake_qchem
from pymatgen.io.qchem_io.inputs import QCInput
import numpy as np

__author__ = "Brandon Wood"
__email__ = "b.wood@berkeley.edu"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")


class TestTorsionPotential(AtomateTest):
    def test_torsion_potential(self):
        # location of test files
        test_tor_files = os.path.join(module_dir, "..", "..", "test_files",
                                      "torsion_wf")
        # define starting molecule and torsion potential workflow object
        initial_qcin = QCInput.from_file(
            os.path.join(test_tor_files, "initial_opt", "mol.qin"))
        initial_mol = initial_qcin.molecule
        atom_indexes = [6, 8, 9, 10]
        angles = [0.0, 90.0, 180.0]
        rem = []
        # add the first rem section
        rem.append({
            "jobtype": "opt",
            "method": "wb97m-v",
            "basis": "def2-tzvppd",
            "gen_scfman": "true",
            "geom_opt_max_cycles": 75,
            "max_scf_cycles": 300,
            "scf_algorithm": "diis",
            "scf_guess": "sad",
            "sym_ignore": "true",
            "symmetry": "false",
            "thresh": 14
        })

        # the second rem section
        rem.append({
            "jobtype": "opt",
            "method": "wb97m-v",
            "basis": "def2-tzvppd",
            "geom_opt_max_cycles": 75,
            "max_scf_cycles": 300,
            "scf_algorithm": "diis",
            "scf_guess": "sad",
            "sym_ignore": "true",
            "symmetry": "false",
            "thresh": 14
        })

        real_wf = get_wf_torsion_potential(
            molecule=initial_mol,
            atom_indexes=atom_indexes,
            angles=angles,
            rem=rem,
            db_file=">>db_file<<")
        # use powerup to replace run with fake run
        # def ref_dirs
        ref_dirs = {
            "initial_opt": os.path.join(test_tor_files, "initial_opt"),
            "opt_0": os.path.join(test_tor_files, "opt_0"),
            "opt_90": os.path.join(test_tor_files, "opt_90"),
            "opt_180": os.path.join(test_tor_files, "opt_180")
        }
        fake_wf = use_fake_qchem(real_wf, ref_dirs)

        self.lp.add_wf(fake_wf)
        rapidfire(
            self.lp,
            fworker=FWorker(env={"db_file": os.path.join(db_dir, "db.json")}))

        wf_test = self.lp.get_wf_by_fw_id(1)
        self.assertTrue(
            all([s == "COMPLETED" for s in wf_test.fw_states.values()]))

        # Checking of the inputs happens in fake_run_qchem so there is no point to retest the inputs
        # Check the output info that gets inserted in the DB
        init_opt = self.get_task_collection().find_one({
            "task_label":
            "initial_opt"
        })
        init_opt_final_mol = Molecule.from_dict(
            init_opt["output"]["optimized_molecule"])
        init_opt_final_e = init_opt["output"]["final_energy"]
        # parse output file
        act_init_opt_out = QCOutput(
            os.path.join(test_tor_files, "initial_opt", "mol.qout"))
        act_init_opt_mol = act_init_opt_out.data[
            "molecule_from_optimized_geometry"]
        act_init_opt_final_e = act_init_opt_out.data["final_energy"]

        np.testing.assert_equal(act_init_opt_mol.species,
                                init_opt_final_mol.species)
        np.testing.assert_allclose(
            act_init_opt_mol.cart_coords,
            init_opt_final_mol.cart_coords,
            atol=0.0001)
        np.testing.assert_equal(act_init_opt_final_e, init_opt_final_e)

        # Optimization of 0 torsion
        opt_0 = self.get_task_collection().find_one({"task_label": "opt_0"})
        opt_0_final_mol = Molecule.from_dict(
            opt_0["output"]["optimized_molecule"])
        opt_0_final_e = opt_0["output"]["final_energy"]
        # parse output file
        act_opt_0_out = QCOutput(
            os.path.join(test_tor_files, "opt_0", "mol.qout"))
        act_opt_0_mol = act_opt_0_out.data["molecule_from_optimized_geometry"]
        act_opt_0_final_e = act_opt_0_out.data["final_energy"]

        np.testing.assert_equal(act_opt_0_mol.species, opt_0_final_mol.species)
        np.testing.assert_allclose(
            act_opt_0_mol.cart_coords,
            opt_0_final_mol.cart_coords,
            atol=0.0001)
        np.testing.assert_equal(act_opt_0_final_e, opt_0_final_e)

        # Optimization of 90 torsion
        opt_90 = self.get_task_collection().find_one({"task_label": "opt_90"})
        opt_90_final_mol = Molecule.from_dict(
            opt_90["output"]["optimized_molecule"])
        opt_90_final_e = opt_90["output"]["final_energy"]
        # parse output file
        act_opt_90_out = QCOutput(
            os.path.join(test_tor_files, "opt_90", "mol.qout"))
        act_opt_90_mol = act_opt_90_out.data[
            "molecule_from_optimized_geometry"]
        act_opt_90_final_e = act_opt_90_out.data["final_energy"]

        np.testing.assert_equal(act_opt_90_mol.species,
                                opt_90_final_mol.species)
        np.testing.assert_allclose(
            act_opt_90_mol.cart_coords,
            opt_90_final_mol.cart_coords,
            atol=0.0001)
        np.testing.assert_equal(act_opt_90_final_e, opt_90_final_e)

        # Optimization of 180 torsion
        opt_180 = self.get_task_collection().find_one({
            "task_label": "opt_180"
        })
        opt_180_final_mol = Molecule.from_dict(
            opt_180["output"]["optimized_molecule"])
        opt_180_final_e = opt_180["output"]["final_energy"]
        # parse output file
        act_opt_180_out = QCOutput(
            os.path.join(test_tor_files, "opt_180", "mol.qout"))
        act_opt_180_mol = act_opt_180_out.data[
            "molecule_from_optimized_geometry"]
        act_opt_180_final_e = act_opt_180_out.data["final_energy"]

        np.testing.assert_equal(act_opt_180_mol.species,
                                opt_180_final_mol.species)
        np.testing.assert_allclose(
            act_opt_180_mol.cart_coords,
            opt_180_final_mol.cart_coords,
            atol=0.0001)
        np.testing.assert_equal(act_opt_180_final_e, opt_180_final_e)


if __name__ == "__main__":
    unittest.main()
