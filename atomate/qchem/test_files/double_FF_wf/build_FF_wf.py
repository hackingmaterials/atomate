# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
import shutil

from atomate.qchem.firetasks.run_calc import RunQChemCustodian
from atomate.qchem.firetasks.write_inputs import WriteInputFromIOSet
from atomate.qchem.firetasks.parse_outputs import QChemToDb
from atomate.qchem.fireworks.core import FrequencyFlatteningOptimizeFW
from atomate.qchem.workflows.base.double_FF_opt import get_wf_double_FF_opt
from fireworks import Firework, Workflow, FWorker
from fireworks.core.launchpad import LaunchPad
from atomate.utils.testing import AtomateTest
from pymatgen.core import Molecule
from pymatgen.io.qchem_io.outputs import QCOutput
from pymatgen.io.qchem_io.inputs import QCInput
import numpy as np

out_file = os.path.join("/global/cscratch1/sd/sblau", "FF_working", "test.qout.opt_0")
qc_out = QCOutput(filename=out_file)
act_mol = qc_out.data["initial_molecule"]
wf = get_wf_double_FF_opt(molecule=act_mol, pcm_dielectric=10.0, max_cores=32, qchem_input_params={"basis_set": "6-311++g**", "overwrite_inputs":{"rem": {"sym_ignore": "true"}}})
lp = LaunchPad.auto_load()
lp.add_wf(wf)
