# coding: utf-8

from __future__ import absolute_import, division, print_function, \
    unicode_literals

# This module defines a workflow for optimizing a molecule first in vacuum and then
# in PCM. Both optimizations will include automatic frequency flattening.

from fireworks import Workflow
from atomate.qchem.fireworks.core import FrequencyFlatteningOptimizeFW
from atomate.utils.utils import get_logger

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "3/19/19"

logger = get_logger(__name__)


def get_wf_vac_2smd_opt(molecule,
                        smd_solvent1,
                        smd_solvent2,
                        linked=True,
                        qchem_input_params=None,
                        name="vac_2smd_opt",
                        suffix="",
                        db_file=">>db_file<<",
                        **kwargs):
    """

    Returns:
        Workflow
    """

    orig_qchem_input_params = qchem_input_params or {}

    # Optimize the molecule in vacuum
    fw1 = FrequencyFlatteningOptimizeFW(
        molecule=molecule,
        name="{}:{}".format(molecule.composition.reduced_formula, "vac_FFopt"+suffix),
        qchem_cmd=">>qchem_cmd<<",
        max_cores=">>max_cores<<",
        qchem_input_params=orig_qchem_input_params,
        linked=linked,
        db_file=db_file)

    smd1_qchem_input_params = {"smd_solvent": smd_solvent1}
    for key in orig_qchem_input_params:
        smd1_qchem_input_params[key] = orig_qchem_input_params[key]

    # Optimize in SMD with solvent1
    fw2 = FrequencyFlatteningOptimizeFW(
        name="{}:{}".format(molecule.composition.reduced_formula, "smd_FFopt"+suffix),
        qchem_cmd=">>qchem_cmd<<",
        max_cores=">>max_cores<<",
        qchem_input_params=smd1_qchem_input_params,
        linked=linked,
        db_file=db_file,
        parents=fw1)

    smd2_qchem_input_params = {"smd_solvent": smd_solvent2}
    for key in orig_qchem_input_params:
        smd2_qchem_input_params[key] = orig_qchem_input_params[key]

    # Optimize in SMD with solvent2
    fw3 = FrequencyFlatteningOptimizeFW(
        name="{}:{}".format(molecule.composition.reduced_formula, "smd_FFopt"+suffix),
        qchem_cmd=">>qchem_cmd<<",
        max_cores=">>max_cores<<",
        qchem_input_params=smd2_qchem_input_params,
        linked=linked,
        db_file=db_file,
        parents=fw1)

    fws = [fw1, fw2, fw3]

    wfname = "{}:{}".format(molecule.composition.reduced_formula, name+suffix)

    return Workflow(fws, name=wfname, **kwargs)
