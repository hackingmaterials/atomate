# coding: utf-8


import os

from atomate.qchem.workflows.base.double_FF_opt import get_wf_double_FF_opt
from fireworks.core.launchpad import LaunchPad
from pymatgen.io.qchem.outputs import QCOutput

out_file = os.path.join("/global/cscratch1/sd/sblau", "FF_working", "test.qout.opt_0")
qc_out = QCOutput(filename=out_file)
act_mol = qc_out.data["initial_molecule"]
wf = get_wf_double_FF_opt(molecule=act_mol, pcm_dielectric=10.0, max_cores=32, qchem_input_params={"basis_set": "6-311++g**", "overwrite_inputs":{"rem": {"sym_ignore": "true"}}})
lp = LaunchPad.auto_load()
lp.add_wf(wf)
