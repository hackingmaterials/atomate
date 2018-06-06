# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
import shutil

from atomate.qchem.firetasks.geo_transformations import RotateTorsion
from atomate.qchem.firetasks.write_inputs import WriteInputFromIOSet
from atomate.qchem.firetasks.parse_outputs import QChemToDb
from atomate.qchem.firetasks.run_calc import RunQChemFake
from fireworks import Firework, Workflow, FWorker
from fireworks.core.rocket_launcher import rapidfire
from atomate.utils.testing import AtomateTest
from pymatgen.core import Molecule
from pymatgen.io.qchem_io.outputs import QCOutput
from pymatgen.io.qchem_io.inputs import QCInput
import numpy as np

__author__ = "Brandon Wood"
__email__ = "b.wood@berkeley.edu"


def use_fake_qchem(original_wf, ref_dirs):
    """
        Replaces all RunQChem commands (i.e. RunQChemDirect, RunQChemCustodian) with RunQChemFake.
        This allows for testing without actually running QChem

        Args:
            original_wf (Workflow)
            ref_dirs (dict): key=firework name, value=path to the reference QChem calculation directory

        Returns:
            Workflow
    """
    for idx_fw, fw in enumerate(original_wf.fws):
        for job_type in ref_dirs.keys():
            if job_type in fw.name:
                for idx_t, t in enumerate(fw.tasks):
                    if "RunQChemCustodian" in str(t) or "RunQChemDirect" in str(t):
                        original_wf.fws[idx_fw].tasks[idx_t] = RunQChemFake(ref_dir=ref_dirs[job_type])

    return original_wf

