import os
import shutil
import unittest

try:
    from unittest.mock import patch
except ImportError:
    from unittest.mock import patch

import numpy as np
from custodian.qchem.handlers import QChemErrorHandler
from custodian.qchem.jobs import QCJob
from pymatgen.io.qchem.outputs import QCOutput

from atomate.qchem.firetasks.run_calc import (
    RunQChemCustodian,
    RunQChemDirect,
    RunQChemFake,
)
from atomate.utils.testing import AtomateTest

__author__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"

module_dir = os.path.dirname(os.path.abspath(__file__))


class TestRunCalcQChem(AtomateTest):
    def setUp(self, lpad=False):
        super().setUp(lpad=False)

    def tearDown(self):
        pass

    def test_RunQChemDirect(self):
        with patch(
            "atomate.qchem.firetasks.run_calc.subprocess.call"
        ) as subprocess_patch:
            with patch("atomate.qchem.firetasks.run_calc.os.putenv") as putenv_patch:
                firetask = RunQChemDirect(
                    qchem_cmd="qchem -slurm -nt 12 co_qc.in mol.qout"
                )
                firetask.run_task(fw_spec={})
                subprocess_patch.assert_called_once()
                putenv_patch.assert_called_once()
                self.assertEqual(
                    subprocess_patch.call_args[0][0],
                    "qchem -slurm -nt 12 co_qc.in mol.qout",
                )
                self.assertEqual(putenv_patch.call_args[0][0], "QCSCRATCH")

    def test_RunQChemCustodian_basic_defaults(self):
        with patch("atomate.qchem.firetasks.run_calc.Custodian") as custodian_patch:
            firetask = RunQChemCustodian(
                qchem_cmd="qchem",
                input_file=os.path.join(
                    module_dir, "..", "..", "test_files", "co_qc.in"
                ),
                max_cores=32,
            )
            firetask.run_task(fw_spec={})
            custodian_patch.assert_called_once()
            self.assertEqual(
                custodian_patch.call_args[0][0][0].as_dict(),
                QChemErrorHandler(
                    input_file=os.path.join(
                        module_dir, "..", "..", "test_files", "co_qc.in"
                    ),
                    output_file="mol.qout",
                ).as_dict(),
            )
            self.assertEqual(
                custodian_patch.call_args[0][1][0].as_dict(),
                QCJob(
                    qchem_command="qchem",
                    max_cores=32,
                    multimode="openmp",
                    input_file=os.path.join(
                        module_dir, "..", "..", "test_files", "co_qc.in"
                    ),
                    output_file="mol.qout",
                ).as_dict(),
            )
            self.assertEqual(
                custodian_patch.call_args[1], {"max_errors": 5, "gzipped_output": True}
            )

    def test_RunQChemCustodian_using_fw_spec_defaults(self):
        with patch("atomate.qchem.firetasks.run_calc.Custodian") as custodian_patch:
            firetask = RunQChemCustodian(
                qchem_cmd=">>qchem_cmd<<",
                input_file=os.path.join(
                    module_dir, "..", "..", "test_files", "co_qc.in"
                ),
                max_cores=">>max_cores<<",
                multimode=">>multimode<<",
            )
            firetask.run_task(
                fw_spec={
                    "_fw_env": {
                        "qchem_cmd": "qchem -slurm",
                        "calc_loc": "/this/is/a/test",
                        "max_cores": 32,
                        "multimode": "openmp",
                        "nboexe": "/path/to/nbo7.i4.exe",
                    }
                }
            )
            custodian_patch.assert_called_once()
            self.assertEqual(
                custodian_patch.call_args[0][0][0].as_dict(),
                QChemErrorHandler(
                    input_file=os.path.join(
                        module_dir, "..", "..", "test_files", "co_qc.in"
                    ),
                    output_file="mol.qout",
                ).as_dict(),
            )
            self.assertEqual(
                custodian_patch.call_args[0][1][0].as_dict(),
                QCJob(
                    qchem_command="qchem -slurm",
                    calc_loc="/this/is/a/test",
                    nboexe="/path/to/nbo7.i4.exe",
                    max_cores=32,
                    multimode="openmp",
                    input_file=os.path.join(
                        module_dir, "..", "..", "test_files", "co_qc.in"
                    ),
                    output_file="mol.qout",
                ).as_dict(),
            )
            self.assertEqual(
                custodian_patch.call_args[1], {"max_errors": 5, "gzipped_output": True}
            )

    def test_RunQChemCustodian_basic_not_defaults(self):
        with patch("atomate.qchem.firetasks.run_calc.Custodian") as custodian_patch:
            firetask = RunQChemCustodian(
                qchem_cmd="qchem -slurm",
                calc_loc="/this/is/a/test",
                multimode="mpi",
                input_file=os.path.join(
                    module_dir, "..", "..", "test_files", "co_qc.in"
                ),
                output_file="this_is_a_test.qout",
                max_cores=4,
                qclog_file="this_is_a_test.qclog",
                suffix="bad_idea",
                save_scratch=True,
                max_errors=137,
                gzipped_output=False,
                handler_group="no_handler",
            )
            firetask.run_task(fw_spec={})
            custodian_patch.assert_called_once()
            self.assertEqual(custodian_patch.call_args[0][0], [])
            self.assertEqual(
                custodian_patch.call_args[0][1][0].as_dict(),
                QCJob(
                    qchem_command="qchem -slurm",
                    calc_loc="/this/is/a/test",
                    multimode="mpi",
                    input_file=os.path.join(
                        module_dir, "..", "..", "test_files", "co_qc.in"
                    ),
                    output_file="this_is_a_test.qout",
                    max_cores=4,
                    qclog_file="this_is_a_test.qclog",
                    suffix="bad_idea",
                    save_scratch=True,
                ).as_dict(),
            )
            self.assertEqual(
                custodian_patch.call_args[1],
                {"max_errors": 137, "gzipped_output": False},
            )

    def test_RunQChemCustodian_using_fw_spec_not_defaults(self):
        with patch("atomate.qchem.firetasks.run_calc.Custodian") as custodian_patch:
            firetask = RunQChemCustodian(
                qchem_cmd=">>qchem_cmd<<",
                multimode=">>multimode<<",
                input_file=os.path.join(
                    module_dir, "..", "..", "test_files", "co_qc.in"
                ),
                output_file="this_is_a_test.qout",
                max_cores=4,
                qclog_file="this_is_a_test.qclog",
                suffix="bad_idea",
                save_scratch=True,
                max_errors=137,
                gzipped_output=False,
                handler_group="no_handler",
            )
            firetask.run_task(
                fw_spec={
                    "_fw_env": {
                        "qchem_cmd": "qchem -slurm",
                        "calc_loc": "/this/is/a/test",
                        "max_cores": 32,
                        "multimode": "mpi",
                    }
                }
            )
            custodian_patch.assert_called_once()
            self.assertEqual(custodian_patch.call_args[0][0], [])
            self.assertEqual(
                custodian_patch.call_args[0][1][0].as_dict(),
                QCJob(
                    qchem_command="qchem -slurm",
                    calc_loc="/this/is/a/test",
                    multimode="mpi",
                    input_file=os.path.join(
                        module_dir, "..", "..", "test_files", "co_qc.in"
                    ),
                    output_file="this_is_a_test.qout",
                    max_cores=4,
                    qclog_file="this_is_a_test.qclog",
                    suffix="bad_idea",
                    save_scratch=True,
                ).as_dict(),
            )
            self.assertEqual(
                custodian_patch.call_args[1],
                {"max_errors": 137, "gzipped_output": False},
            )

    def test_RunQChemCustodian_FF_basic_defaults(self):
        with patch("atomate.qchem.firetasks.run_calc.Custodian") as custodian_patch:
            with patch(
                "atomate.qchem.firetasks.run_calc.QCJob.opt_with_frequency_flattener"
            ) as FF_patch:
                firetask = RunQChemCustodian(
                    qchem_cmd="qchem",
                    max_cores=32,
                    input_file=os.path.join(
                        module_dir,
                        "..",
                        "..",
                        "test_files",
                        "FF_before_run",
                        "test.qin",
                    ),
                    output_file=os.path.join(
                        module_dir,
                        "..",
                        "..",
                        "test_files",
                        "FF_before_run",
                        "test.qout",
                    ),
                    job_type="opt_with_frequency_flattener",
                )
                firetask.run_task(fw_spec={})
                custodian_patch.assert_called_once()
                self.assertEqual(
                    FF_patch.call_args[1],
                    {
                        "qchem_command": "qchem",
                        "multimode": "openmp",
                        "input_file": os.path.join(
                            module_dir,
                            "..",
                            "..",
                            "test_files",
                            "FF_before_run",
                            "test.qin",
                        ),
                        "output_file": os.path.join(
                            module_dir,
                            "..",
                            "..",
                            "test_files",
                            "FF_before_run",
                            "test.qout",
                        ),
                        "qclog_file": "mol.qclog",
                        "max_iterations": 10,
                        "linked": True,
                        "save_final_scratch": False,
                        "max_cores": 32,
                        "calc_loc": None,
                        "freq_before_opt": False,
                        "transition_state": False,
                        "nboexe": None,
                    },
                )

    def test_RunQChemCustodian_FF_using_fw_spec_defaults(self):
        with patch("atomate.qchem.firetasks.run_calc.Custodian") as custodian_patch:
            with patch(
                "atomate.qchem.firetasks.run_calc.QCJob.opt_with_frequency_flattener"
            ) as FF_patch:
                firetask = RunQChemCustodian(
                    qchem_cmd=">>qchem_cmd<<",
                    max_cores=">>max_cores<<",
                    multimode=">>multimode<<",
                    input_file=os.path.join(
                        module_dir,
                        "..",
                        "..",
                        "test_files",
                        "FF_before_run",
                        "test.qin",
                    ),
                    output_file=os.path.join(
                        module_dir,
                        "..",
                        "..",
                        "test_files",
                        "FF_before_run",
                        "test.qout",
                    ),
                    job_type="opt_with_frequency_flattener",
                )
                firetask.run_task(
                    fw_spec={
                        "_fw_env": {
                            "qchem_cmd": "qchem -slurm",
                            "calc_loc": "/this/is/a/test",
                            "max_cores": 32,
                            "multimode": "openmp",
                        }
                    }
                )
                custodian_patch.assert_called_once()
                self.assertEqual(
                    FF_patch.call_args[1],
                    {
                        "qchem_command": "qchem -slurm",
                        "multimode": "openmp",
                        "input_file": os.path.join(
                            module_dir,
                            "..",
                            "..",
                            "test_files",
                            "FF_before_run",
                            "test.qin",
                        ),
                        "output_file": os.path.join(
                            module_dir,
                            "..",
                            "..",
                            "test_files",
                            "FF_before_run",
                            "test.qout",
                        ),
                        "qclog_file": "mol.qclog",
                        "max_iterations": 10,
                        "linked": True,
                        "calc_loc": "/this/is/a/test",
                        "save_final_scratch": False,
                        "max_cores": 32,
                        "freq_before_opt": False,
                        "transition_state": False,
                        "nboexe": None,
                    },
                )

    def test_RunQChemCustodian_FF_basic_not_defaults(self):
        with patch("atomate.qchem.firetasks.run_calc.Custodian") as custodian_patch:
            with patch(
                "atomate.qchem.firetasks.run_calc.QCJob.opt_with_frequency_flattener"
            ) as FF_patch:
                firetask = RunQChemCustodian(
                    qchem_cmd="qchem -slurm",
                    calc_loc="/this/is/a/test",
                    nboexe="/path/to/nbo7.i4.exe",
                    input_file=os.path.join(
                        module_dir,
                        "..",
                        "..",
                        "test_files",
                        "FF_before_run",
                        "test.qin",
                    ),
                    output_file=os.path.join(
                        module_dir,
                        "..",
                        "..",
                        "test_files",
                        "FF_before_run",
                        "test.qout",
                    ),
                    job_type="opt_with_frequency_flattener",
                    max_cores=4,
                    qclog_file="this_is_a_test.qclog",
                    suffix="bad_idea",
                    save_scratch=True,
                    max_errors=137,
                    linked=False,
                    gzipped_output=False,
                    handler_group="no_handler",
                    max_iterations=1029,
                    max_molecule_perturb_scale=0.5,
                    multimode="mpi",
                    transition_state=True,
                    freq_before_opt=True,
                )
                firetask.run_task(fw_spec={})
                custodian_patch.assert_called_once()
                self.assertEqual(custodian_patch.call_args[0][0], [])
                self.assertEqual(
                    FF_patch.call_args[1],
                    {
                        "qchem_command": "qchem -slurm",
                        "multimode": "mpi",
                        "input_file": os.path.join(
                            module_dir,
                            "..",
                            "..",
                            "test_files",
                            "FF_before_run",
                            "test.qin",
                        ),
                        "output_file": os.path.join(
                            module_dir,
                            "..",
                            "..",
                            "test_files",
                            "FF_before_run",
                            "test.qout",
                        ),
                        "qclog_file": "this_is_a_test.qclog",
                        "max_iterations": 1029,
                        "max_molecule_perturb_scale": 0.5,
                        "linked": False,
                        "calc_loc": "/this/is/a/test",
                        "nboexe": "/path/to/nbo7.i4.exe",
                        "save_final_scratch": True,
                        "max_cores": 4,
                        "freq_before_opt": True,
                        "transition_state": True,
                    },
                )

    def test_RunQChemCustodian_FF_using_fw_spec_not_defaults(self):
        with patch("atomate.qchem.firetasks.run_calc.Custodian") as custodian_patch:
            with patch(
                "atomate.qchem.firetasks.run_calc.QCJob.opt_with_frequency_flattener"
            ) as FF_patch:
                firetask = RunQChemCustodian(
                    qchem_cmd=">>qchem_cmd<<",
                    input_file=os.path.join(
                        module_dir,
                        "..",
                        "..",
                        "test_files",
                        "FF_before_run",
                        "test.qin",
                    ),
                    output_file=os.path.join(
                        module_dir,
                        "..",
                        "..",
                        "test_files",
                        "FF_before_run",
                        "test.qout",
                    ),
                    job_type="opt_with_frequency_flattener",
                    max_cores=4,
                    qclog_file="this_is_a_test.qclog",
                    suffix="bad_idea",
                    save_scratch=True,
                    max_errors=137,
                    gzipped_output=False,
                    handler_group="no_handler",
                    linked=False,
                    max_iterations=1029,
                    max_molecule_perturb_scale=0.5,
                    multimode=">>multimode<<",
                    transition_state=True,
                    freq_before_opt=True,
                )
                firetask.run_task(
                    fw_spec={
                        "_fw_env": {
                            "qchem_cmd": "qchem -slurm",
                            "calc_loc": "/this/is/a/test",
                            "max_cores": 32,
                            "multimode": "mpi",
                        }
                    }
                )
                custodian_patch.assert_called_once()
                self.assertEqual(custodian_patch.call_args[0][0], [])
                self.assertEqual(
                    FF_patch.call_args[1],
                    {
                        "qchem_command": "qchem -slurm",
                        "multimode": "mpi",
                        "input_file": os.path.join(
                            module_dir,
                            "..",
                            "..",
                            "test_files",
                            "FF_before_run",
                            "test.qin",
                        ),
                        "output_file": os.path.join(
                            module_dir,
                            "..",
                            "..",
                            "test_files",
                            "FF_before_run",
                            "test.qout",
                        ),
                        "qclog_file": "this_is_a_test.qclog",
                        "max_iterations": 1029,
                        "max_molecule_perturb_scale": 0.5,
                        "linked": False,
                        "calc_loc": "/this/is/a/test",
                        "nboexe": None,
                        "save_final_scratch": True,
                        "max_cores": 4,
                        "freq_before_opt": True,
                        "transition_state": True,
                    },
                )


class TestFakeRunQChem(AtomateTest):
    def setUp(self, lpad=False):
        os.makedirs(os.path.join(module_dir, "..", "fake_run"))
        shutil.copyfile(
            os.path.join(module_dir, "..", "..", "test_files", "real_run", "mol.qin"),
            os.path.join(module_dir, "..", "fake_run", "mol.qin"),
        )
        super().setUp(lpad=False)

    def tearDown(self):
        shutil.rmtree(os.path.join(module_dir, "..", "fake_run"))

    def test_RunQChemFake(self):
        os.chdir(os.path.join(module_dir, "..", "fake_run"))
        firetask = RunQChemFake(
            ref_dir=os.path.join(module_dir, "..", "..", "test_files", "real_run")
        )
        firetask.run_task(fw_spec={})
        ref_out = QCOutput(
            os.path.join(module_dir, "..", "..", "test_files", "real_run", "mol.qout")
        ).data
        this_out = QCOutput("mol.qout").data
        for key in ref_out:
            if key == "dipoles":
                self.assertEqual(ref_out[key]["total"], this_out[key]["total"])
                np.testing.assert_array_equal(ref_out[key]["dipole"], this_out[key]["dipole"])
            else:
                try:
                    self.assertEqual(ref_out[key], this_out[key])
                except ValueError:
                    np.testing.assert_array_equal(ref_out[key], this_out[key])
        for filename in os.listdir(
            os.path.join(module_dir, "..", "..", "test_files", "real_run")
        ):
            self.assertEqual(os.path.isfile(filename), True)
        os.chdir(module_dir)


if __name__ == "__main__":
    unittest.main()
