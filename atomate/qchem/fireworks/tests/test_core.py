import os
import unittest
from itertools import chain

import numpy as np
from pymatgen.io.qchem.outputs import QCOutput
from pymatgen.core.structure import Molecule
from pymatgen.core.sites import Site

from atomate.qchem.firetasks.critic2 import ProcessCritic2, RunCritic2
from atomate.qchem.firetasks.fragmenter import FragmentMolecule
from atomate.qchem.firetasks.geo_transformations import PerturbGeometry
from atomate.qchem.firetasks.parse_outputs import QChemToDb
from atomate.qchem.firetasks.parse_outputs import ProtCalcToDb
from atomate.qchem.firetasks.run_calc import RunQChemCustodian
from atomate.qchem.firetasks.write_inputs import WriteInputFromIOSet
from atomate.qchem.fireworks.core import (
    CubeAndCritic2FW,
    FragmentFW,
    FrequencyFlatteningOptimizeFW,
    FrequencyFlatteningTransitionStateFW,
    FrequencyFW,
    OptimizeFW,
    PESScanFW,
    ProtonEnergyFW,
    SinglePointFW,
    TransitionStateFW,
)
from atomate.utils.testing import AtomateTest

__author__ = "Samuel Blau, Evan Spotte-Smith"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "5/23/18"
__credits__ = "Brandon Wood, Shyam Dwaraknath"

module_dir = os.path.dirname(os.path.abspath(__file__))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
db_file = os.path.join(db_dir, "db.json")


class TestCore(AtomateTest):
    def setUp(self, lpad=False):
        out_file = os.path.join(
            module_dir, "..", "..", "test_files", "FF_working", "test.qout.opt_0"
        )
        qc_out = QCOutput(filename=out_file)
        self.act_mol = qc_out.data["initial_molecule"]

        self.maxDiff = None

        super().setUp(lpad=False)

    def tearDown(self):
        pass

    def test_SinglePointFW_defaults(self):
        firework = SinglePointFW(molecule=self.act_mol)
        self.assertEqual(
            firework.tasks[0].as_dict(),
            WriteInputFromIOSet(
                molecule=self.act_mol,
                qchem_input_set="SinglePointSet",
                input_file="mol.qin",
                qchem_input_params={},
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[1].as_dict(),
            RunQChemCustodian(
                qchem_cmd=">>qchem_cmd<<",
                multimode=">>multimode<<",
                input_file="mol.qin",
                output_file="mol.qout",
                max_cores=">>max_cores<<",
                max_errors=5,
                job_type="normal",
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[2].as_dict(),
            QChemToDb(
                db_file=None,
                input_file="mol.qin",
                output_file="mol.qout",
                additional_fields={"task_label": "single point"},
            ).as_dict(),
        )
        self.assertEqual(firework.parents, [])
        self.assertEqual(firework.name, "single point")

    def test_SinglePointFW_not_defaults(self):
        firework = SinglePointFW(
            molecule=self.act_mol,
            name="special single point",
            qchem_cmd="qchem -slurm",
            multimode="mpi",
            max_cores=12,
            qchem_input_params={"pcm_dielectric": 10.0},
            db_file=db_file,
            parents=None,
        )
        self.assertEqual(
            firework.tasks[0].as_dict(),
            WriteInputFromIOSet(
                molecule=self.act_mol,
                qchem_input_set="SinglePointSet",
                input_file="mol.qin",
                qchem_input_params={"pcm_dielectric": 10.0},
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[1].as_dict(),
            RunQChemCustodian(
                qchem_cmd="qchem -slurm",
                multimode="mpi",
                input_file="mol.qin",
                output_file="mol.qout",
                max_cores=12,
                max_errors=5,
                job_type="normal",
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[2].as_dict(),
            QChemToDb(
                db_file=db_file,
                input_file="mol.qin",
                output_file="mol.qout",
                additional_fields={"task_label": "special single point"},
            ).as_dict(),
        )
        self.assertEqual(firework.parents, [])
        self.assertEqual(firework.name, "special single point")

    def test_ProtonEnergyFW(self):

        H_site_1_H2O = Site("H", [0.18338, 2.20176, 0.01351])
        H_site_2_H2O = Site("H", [-1.09531, 1.61602, 0.70231])
        O_site_H2O = Site("O", [-0.80595, 2.22952, -0.01914])
        H2O_molecule = Molecule.from_sites([H_site_1_H2O, H_site_2_H2O, O_site_H2O])

        H_site_1_H3O = Site("H", [0.11550, 2.34733, 0.00157])
        H_site_2_H3O = Site("H", [-1.17463, 1.77063, 0.67652])
        H_site_3_H3O = Site("H", [-1.29839, 2.78012, -0.51436])
        O_site_H3O = Site("O", [-0.78481, 1.99137, -0.20661])
        H3O_ion = Molecule.from_sites(
            [H_site_1_H3O, H_site_2_H3O, H_site_3_H3O, O_site_H3O]
        )

        H2O_molecule.set_charge_and_spin(0, 1)
        H3O_ion.set_charge_and_spin(1, 1)

        firework = ProtonEnergyFW(qchem_input_params={"smd_solvent": "water"})
        self.assertEqual(
            firework.tasks[0].as_dict(),
            WriteInputFromIOSet(
                molecule=H2O_molecule,
                qchem_input_set="OptSet",
                input_file="water.qin",
                qchem_input_params={"smd_solvent": "water"},
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[1].as_dict(),
            RunQChemCustodian(
                qchem_cmd=">>qchem_cmd<<",
                multimode=">>multimode<<",
                input_file="water.qin",
                output_file="water.qout",
                max_cores=">>max_cores<<",
                max_errors=5,
                job_type="normal",
                gzipped_output=False,
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[2].as_dict(),
            WriteInputFromIOSet(
                molecule=H3O_ion,
                qchem_input_set="OptSet",
                input_file="hydronium.qin",
                qchem_input_params={"smd_solvent": "water"},
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[3].as_dict(),
            RunQChemCustodian(
                qchem_cmd=">>qchem_cmd<<",
                multimode=">>multimode<<",
                input_file="hydronium.qin",
                output_file="hydronium.qout",
                max_cores=">>max_cores<<",
                max_errors=5,
                job_type="normal",
                gzipped_output=False,
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[4].as_dict(),
            ProtCalcToDb(
                db_file=None,
                input_file_H2O="water.qin",
                output_file_H2O="water.qout",
                input_file_H3O="hydronium.qin",
                output_file_H3O="hydronium.qout",
                additional_fields={"task_label": "proton electronic energy"},
            ).as_dict(),
        )
        self.assertEqual(firework.parents, [])
        self.assertEqual(firework.name, "proton electronic energy")

    def test_OptimizeFW_defaults(self):
        firework = OptimizeFW(molecule=self.act_mol)
        self.assertEqual(
            firework.tasks[0].as_dict(),
            WriteInputFromIOSet(
                molecule=self.act_mol,
                qchem_input_set="OptSet",
                input_file="mol.qin",
                qchem_input_params={},
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[1].as_dict(),
            RunQChemCustodian(
                qchem_cmd=">>qchem_cmd<<",
                multimode=">>multimode<<",
                input_file="mol.qin",
                output_file="mol.qout",
                max_cores=">>max_cores<<",
                max_errors=20,
                job_type="normal",
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[2].as_dict(),
            QChemToDb(
                db_file=None,
                input_file="mol.qin",
                output_file="mol.qout",
                additional_fields={"task_label": "structure optimization"},
            ).as_dict(),
        )
        self.assertEqual(firework.parents, [])
        self.assertEqual(firework.name, "structure optimization")

    def test_OptimizeFW_not_defaults(self):
        firework = OptimizeFW(
            molecule=self.act_mol,
            name="special structure optimization",
            qchem_cmd="qchem -slurm",
            multimode="mpi",
            max_cores=12,
            qchem_input_params={"pcm_dielectric": 10.0},
            db_file=db_file,
            parents=None,
        )
        self.assertEqual(
            firework.tasks[0].as_dict(),
            WriteInputFromIOSet(
                molecule=self.act_mol,
                qchem_input_set="OptSet",
                input_file="mol.qin",
                qchem_input_params={"pcm_dielectric": 10.0},
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[1].as_dict(),
            RunQChemCustodian(
                qchem_cmd="qchem -slurm",
                multimode="mpi",
                input_file="mol.qin",
                output_file="mol.qout",
                max_cores=12,
                max_errors=20,
                job_type="normal",
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[2].as_dict(),
            QChemToDb(
                db_file=db_file,
                input_file="mol.qin",
                output_file="mol.qout",
                additional_fields={"task_label": "special structure optimization"},
            ).as_dict(),
        )
        self.assertEqual(firework.parents, [])
        self.assertEqual(firework.name, "special structure optimization")

    def test_TransitionStateFW_defaults(self):
        firework = TransitionStateFW(molecule=self.act_mol)
        self.assertEqual(
            firework.tasks[0].as_dict(),
            WriteInputFromIOSet(
                molecule=self.act_mol,
                qchem_input_set="TransitionStateSet",
                input_file="mol.qin",
                qchem_input_params={},
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[1].as_dict(),
            RunQChemCustodian(
                qchem_cmd=">>qchem_cmd<<",
                multimode=">>multimode<<",
                input_file="mol.qin",
                output_file="mol.qout",
                max_cores=">>max_cores<<",
                max_errors=5,
                job_type="normal",
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[2].as_dict(),
            QChemToDb(
                db_file=None,
                input_file="mol.qin",
                output_file="mol.qout",
                additional_fields={
                    "task_label": "transition state structure optimization"
                },
            ).as_dict(),
        )
        self.assertEqual(firework.parents, [])
        self.assertEqual(firework.name, "transition state structure optimization")

    def test_TransitionStateFW_not_defaults(self):
        firework = TransitionStateFW(
            molecule=self.act_mol,
            name="special transition state structure optimization",
            qchem_cmd="qchem -slurm",
            multimode="mpi",
            max_cores=12,
            qchem_input_params={"pcm_dielectric": 10.0},
            db_file=db_file,
            parents=None,
        )
        self.assertEqual(
            firework.tasks[0].as_dict(),
            WriteInputFromIOSet(
                molecule=self.act_mol,
                qchem_input_set="TransitionStateSet",
                input_file="mol.qin",
                qchem_input_params={"pcm_dielectric": 10.0},
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[1].as_dict(),
            RunQChemCustodian(
                qchem_cmd="qchem -slurm",
                multimode="mpi",
                input_file="mol.qin",
                output_file="mol.qout",
                max_cores=12,
                max_errors=5,
                job_type="normal",
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[2].as_dict(),
            QChemToDb(
                db_file=db_file,
                input_file="mol.qin",
                output_file="mol.qout",
                additional_fields={
                    "task_label": "special transition state structure optimization"
                },
            ).as_dict(),
        )
        self.assertEqual(firework.parents, [])
        self.assertEqual(
            firework.name, "special transition state structure optimization"
        )

    def test_FrequencyFW_defaults(self):
        firework = FrequencyFW(molecule=self.act_mol)
        self.assertEqual(
            firework.tasks[0].as_dict(),
            WriteInputFromIOSet(
                molecule=self.act_mol,
                qchem_input_set="FreqSet",
                input_file="mol.qin",
                qchem_input_params={},
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[1].as_dict(),
            RunQChemCustodian(
                qchem_cmd=">>qchem_cmd<<",
                multimode=">>multimode<<",
                input_file="mol.qin",
                output_file="mol.qout",
                max_cores=">>max_cores<<",
                max_errors=5,
                job_type="normal",
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[2].as_dict(),
            QChemToDb(
                db_file=None,
                input_file="mol.qin",
                output_file="mol.qout",
                additional_fields={"task_label": "frequency calculation"},
            ).as_dict(),
        )
        self.assertEqual(firework.parents, [])
        self.assertEqual(firework.name, "frequency calculation")

    def test_FrequencyFW_not_defaults(self):
        firework = FrequencyFW(
            molecule=self.act_mol,
            name="special frequency analysis",
            qchem_cmd="qchem -slurm",
            multimode="mpi",
            max_cores=12,
            qchem_input_params={"pcm_dielectric": 10.0},
            db_file=db_file,
            parents=None,
        )
        self.assertEqual(
            firework.tasks[0].as_dict(),
            WriteInputFromIOSet(
                molecule=self.act_mol,
                qchem_input_set="FreqSet",
                input_file="mol.qin",
                qchem_input_params={"pcm_dielectric": 10.0},
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[1].as_dict(),
            RunQChemCustodian(
                qchem_cmd="qchem -slurm",
                multimode="mpi",
                input_file="mol.qin",
                output_file="mol.qout",
                max_cores=12,
                max_errors=5,
                job_type="normal",
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[2].as_dict(),
            QChemToDb(
                db_file=db_file,
                input_file="mol.qin",
                output_file="mol.qout",
                additional_fields={"task_label": "special frequency analysis"},
            ).as_dict(),
        )
        self.assertEqual(firework.parents, [])
        self.assertEqual(firework.name, "special frequency analysis")

    def test_FrequencyFlatteningOptimizeFW_defaults(self):
        firework = FrequencyFlatteningOptimizeFW(molecule=self.act_mol)
        self.assertEqual(
            firework.tasks[0].as_dict(),
            WriteInputFromIOSet(
                molecule=self.act_mol,
                qchem_input_set="OptSet",
                input_file="mol.qin",
                qchem_input_params={},
                prev_hess=None,
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[1].as_dict(),
            RunQChemCustodian(
                qchem_cmd=">>qchem_cmd<<",
                multimode=">>multimode<<",
                input_file="mol.qin",
                output_file="mol.qout",
                max_cores=">>max_cores<<",
                job_type="opt_with_frequency_flattener",
                max_iterations=10,
                max_molecule_perturb_scale=0.3,
                linked=True,
                freq_before_opt=False,
                max_errors=20,
                save_scratch=True,
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[2].as_dict(),
            QChemToDb(
                db_file=None,
                input_file="mol.qin",
                output_file="mol.qout",
                parse_hess_file=True,
                additional_fields={
                    "task_label": "frequency flattening structure optimization",
                    "special_run_type": "frequency_flattener",
                    "linked": True,
                },
            ).as_dict(),
        )
        self.assertEqual(firework.parents, [])
        self.assertEqual(firework.name, "frequency flattening structure optimization")

    def test_FrequencyFlatteningOptimizeFW_not_defaults(self):
        firework = FrequencyFlatteningOptimizeFW(
            molecule=self.act_mol,
            name="special frequency flattening structure optimization",
            qchem_cmd="qchem -slurm",
            multimode="mpi",
            max_cores=12,
            qchem_input_params={"pcm_dielectric": 10.0},
            max_iterations=5,
            max_molecule_perturb_scale=0.2,
            linked=False,
            freq_before_opt=True,
            db_file=db_file,
            perturb_geometry=True,
            mode=np.zeros((len(self.act_mol), 3)),
            scale=0.2,
            parents=None,
            max_errors=20,
        )
        self.assertEqual(
            firework.tasks[0].as_dict(),
            PerturbGeometry(
                molecule=self.act_mol, mode=np.zeros((len(self.act_mol), 3)), scale=0.2
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[1].as_dict(),
            WriteInputFromIOSet(
                molecule=None,
                qchem_input_set="FreqSet",
                input_file="mol.qin",
                qchem_input_params={"pcm_dielectric": 10.0},
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[2].as_dict(),
            RunQChemCustodian(
                qchem_cmd="qchem -slurm",
                multimode="mpi",
                input_file="mol.qin",
                output_file="mol.qout",
                max_cores=12,
                job_type="opt_with_frequency_flattener",
                max_iterations=5,
                max_molecule_perturb_scale=0.2,
                linked=False,
                freq_before_opt=True,
                max_errors=20,
                save_scratch=True,
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[3].as_dict(),
            QChemToDb(
                db_file=db_file,
                input_file="mol.qin",
                output_file="mol.qout",
                parse_hess_file=True,
                additional_fields={
                    "task_label": "special frequency flattening structure optimization",
                    "special_run_type": "frequency_flattener",
                    "linked": False,
                },
            ).as_dict(),
        )
        self.assertEqual(firework.parents, [])
        self.assertEqual(
            firework.name, "special frequency flattening structure optimization"
        )

    def test_FrequencyFlatteningTransitionStateFW_defaults(self):
        firework = FrequencyFlatteningTransitionStateFW(molecule=self.act_mol)
        self.assertEqual(
            firework.tasks[0].as_dict(),
            WriteInputFromIOSet(
                molecule=self.act_mol,
                qchem_input_set="FreqSet",
                input_file="mol.qin",
                qchem_input_params={},
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[1].as_dict(),
            RunQChemCustodian(
                qchem_cmd=">>qchem_cmd<<",
                multimode=">>multimode<<",
                input_file="mol.qin",
                output_file="mol.qout",
                max_cores=">>max_cores<<",
                job_type="opt_with_frequency_flattener",
                max_iterations=3,
                max_molecule_perturb_scale=0.3,
                transition_state=True,
                freq_before_opt=True,
                linked=True,
                max_errors=5,
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[2].as_dict(),
            QChemToDb(
                db_file=None,
                input_file="mol.qin",
                output_file="mol.qout",
                additional_fields={
                    "task_label": "frequency flattening transition state optimization",
                    "special_run_type": "ts_frequency_flattener",
                    "linked": True,
                },
                runs=["freq_pre"]
                + list(
                    chain.from_iterable(
                        [["ts_" + str(ii), "freq_" + str(ii)] for ii in range(10)]
                    )
                ),
            ).as_dict(),
        )
        self.assertEqual(firework.parents, [])
        self.assertEqual(
            firework.name, "frequency flattening transition state optimization"
        )

    def test_FrequencyFlatteningTransitionStateFW_not_defaults(self):
        self.maxDiff = None
        firework = FrequencyFlatteningTransitionStateFW(
            molecule=self.act_mol,
            name="special frequency flattening transition state optimization",
            qchem_cmd="qchem -slurm",
            multimode="mpi",
            max_cores=12,
            qchem_input_params={"pcm_dielectric": 10.0},
            max_iterations=5,
            max_molecule_perturb_scale=0.2,
            linked=False,
            freq_before_opt=False,
            perturb_geometry=True,
            mode=np.zeros((len(self.act_mol), 3)),
            scale=0.2,
            db_file=db_file,
            parents=None,
            max_errors=5,
        )
        self.assertEqual(
            firework.tasks[0].as_dict(),
            PerturbGeometry(
                molecule=self.act_mol, mode=np.zeros((len(self.act_mol), 3)), scale=0.2
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[1].as_dict(),
            WriteInputFromIOSet(
                molecule=None,
                qchem_input_set="TransitionStateSet",
                input_file="mol.qin",
                qchem_input_params={"pcm_dielectric": 10.0},
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[2].as_dict(),
            RunQChemCustodian(
                qchem_cmd="qchem -slurm",
                multimode="mpi",
                input_file="mol.qin",
                output_file="mol.qout",
                max_cores=12,
                job_type="opt_with_frequency_flattener",
                max_iterations=5,
                max_molecule_perturb_scale=0.2,
                transition_state=True,
                linked=False,
                freq_before_opt=False,
                max_errors=5,
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[3].as_dict(),
            QChemToDb(
                db_file=db_file,
                input_file="mol.qin",
                output_file="mol.qout",
                additional_fields={
                    "task_label": "special frequency flattening transition state optimization",
                    "special_run_type": "ts_frequency_flattener",
                    "linked": False,
                },
                runs=list(
                    chain.from_iterable(
                        [["ts_" + str(ii), "freq_" + str(ii)] for ii in range(10)]
                    )
                ),
            ).as_dict(),
        )
        self.assertEqual(firework.parents, [])
        self.assertEqual(
            firework.name, "special frequency flattening transition state optimization"
        )

    def test_PESScanFW_defaults(self):
        firework = PESScanFW(
            molecule=self.act_mol, scan_variables={"stre": ["0 1 1.5 2.0 0.01"]}
        )
        self.assertEqual(
            firework.tasks[0].as_dict(),
            WriteInputFromIOSet(
                molecule=self.act_mol,
                qchem_input_set="PESScanSet",
                input_file="mol.qin",
                qchem_input_params={"scan_variables": {"stre": ["0 1 1.5 2.0 0.01"]}},
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[1].as_dict(),
            RunQChemCustodian(
                qchem_cmd=">>qchem_cmd<<",
                multimode=">>multimode<<",
                input_file="mol.qin",
                output_file="mol.qout",
                max_cores=">>max_cores<<",
                max_errors=5,
                job_type="normal",
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[2].as_dict(),
            QChemToDb(
                db_file=None,
                input_file="mol.qin",
                output_file="mol.qout",
                additional_fields={"task_label": "potential energy surface scan"},
            ).as_dict(),
        )
        self.assertEqual(firework.parents, [])
        self.assertEqual(firework.name, "potential energy surface scan")

    def test_PESScanFW_not_defaults(self):
        firework = PESScanFW(
            molecule=self.act_mol,
            name="special potential energy surface scan",
            qchem_cmd="qchem -slurm",
            multimode="mpi",
            max_cores=12,
            qchem_input_params={"pcm_dielectric": 10.0},
            db_file=db_file,
            parents=None,
            scan_variables={"stre": ["0 1 1.5 2.0 0.01"]},
        )
        self.assertEqual(
            firework.tasks[0].as_dict(),
            WriteInputFromIOSet(
                molecule=self.act_mol,
                qchem_input_set="PESScanSet",
                input_file="mol.qin",
                qchem_input_params={
                    "pcm_dielectric": 10.0,
                    "scan_variables": {"stre": ["0 1 1.5 2.0 0.01"]},
                },
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[1].as_dict(),
            RunQChemCustodian(
                qchem_cmd="qchem -slurm",
                multimode="mpi",
                input_file="mol.qin",
                output_file="mol.qout",
                max_cores=12,
                job_type="normal",
                max_errors=5,
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[2].as_dict(),
            QChemToDb(
                db_file=db_file,
                input_file="mol.qin",
                output_file="mol.qout",
                additional_fields={
                    "task_label": "special potential energy surface scan"
                },
            ).as_dict(),
        )
        self.assertEqual(firework.parents, [])
        self.assertEqual(firework.name, "special potential energy surface scan")

    def test_FragmentFW_defaults(self):
        firework = FragmentFW(molecule=self.act_mol)
        self.assertEqual(
            firework.tasks[0].as_dict(),
            FragmentMolecule(
                molecule=self.act_mol,
                depth=1,
                open_rings=True,
                additional_charges=[],
                do_triplets=True,
                linked=False,
                qchem_input_params={},
                db_file=None,
                check_db=True,
            ).as_dict(),
        )
        self.assertEqual(firework.parents, [])
        self.assertEqual(firework.name, "fragment and optimize")

    def test_FragmentFW_not_defaults(self):
        firework = FragmentFW(
            molecule=self.act_mol,
            depth=0,
            open_rings=False,
            additional_charges=[2],
            do_triplets=False,
            linked=True,
            name="fragmenting a thing",
            qchem_input_params={"pcm_dielectric": 10.0},
            db_file=db_file,
            check_db=False,
        )
        self.assertEqual(
            firework.tasks[0].as_dict(),
            FragmentMolecule(
                molecule=self.act_mol,
                depth=0,
                open_rings=False,
                additional_charges=[2],
                do_triplets=False,
                linked=True,
                qchem_input_params={"pcm_dielectric": 10.0},
                db_file=db_file,
                check_db=False,
            ).as_dict(),
        )
        self.assertEqual(firework.parents, [])
        self.assertEqual(firework.name, "fragmenting a thing")

    def test_CubeAndCritic2FW_defaults(self):
        firework = CubeAndCritic2FW(molecule=self.act_mol)
        self.assertEqual(
            firework.tasks[0].as_dict(),
            WriteInputFromIOSet(
                molecule=self.act_mol,
                qchem_input_set="SinglePointSet",
                input_file="mol.qin",
                qchem_input_params={"plot_cubes": True},
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[1].as_dict(),
            RunQChemCustodian(
                qchem_cmd=">>qchem_cmd<<",
                multimode=">>multimode<<",
                input_file="mol.qin",
                output_file="mol.qout",
                max_cores=">>max_cores<<",
                max_errors=5,
                job_type="normal",
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[2].as_dict(),
            RunCritic2(molecule=self.act_mol, cube_file="dens.0.cube.gz").as_dict(),
        )
        self.assertEqual(
            firework.tasks[3].as_dict(), ProcessCritic2(molecule=self.act_mol).as_dict()
        )
        self.assertEqual(
            firework.tasks[4].as_dict(),
            QChemToDb(
                db_file=None,
                input_file="mol.qin",
                output_file="mol.qout",
                additional_fields={"task_label": "cube and critic2"},
            ).as_dict(),
        )
        self.assertEqual(firework.parents, [])
        self.assertEqual(firework.name, "cube and critic2")

    def test_CubeAndCritic2FW_not_defaults(self):
        firework = CubeAndCritic2FW(
            molecule=self.act_mol,
            name="special cube and critic2",
            qchem_cmd="qchem -slurm",
            multimode="mpi",
            max_cores=12,
            qchem_input_params={"pcm_dielectric": 10.0},
            db_file=db_file,
            parents=None,
        )
        self.assertEqual(
            firework.tasks[0].as_dict(),
            WriteInputFromIOSet(
                molecule=self.act_mol,
                qchem_input_set="SinglePointSet",
                input_file="mol.qin",
                qchem_input_params={"pcm_dielectric": 10.0, "plot_cubes": True},
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[1].as_dict(),
            RunQChemCustodian(
                qchem_cmd="qchem -slurm",
                multimode="mpi",
                input_file="mol.qin",
                output_file="mol.qout",
                max_cores=12,
                max_errors=5,
                job_type="normal",
            ).as_dict(),
        )
        self.assertEqual(
            firework.tasks[2].as_dict(),
            RunCritic2(molecule=self.act_mol, cube_file="dens.0.cube.gz").as_dict(),
        )
        self.assertEqual(
            firework.tasks[3].as_dict(), ProcessCritic2(molecule=self.act_mol).as_dict()
        )
        self.assertEqual(
            firework.tasks[4].as_dict(),
            QChemToDb(
                db_file=db_file,
                input_file="mol.qin",
                output_file="mol.qout",
                additional_fields={"task_label": "special cube and critic2"},
            ).as_dict(),
        )
        self.assertEqual(firework.parents, [])
        self.assertEqual(firework.name, "special cube and critic2")


if __name__ == "__main__":
    unittest.main()
