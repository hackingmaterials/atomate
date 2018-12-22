# coding: utf-8

from __future__ import absolute_import, division, print_function, \
    unicode_literals

# This module defines a workflow for optimizing a molecule first in vacuum and then
# in PCM. Both optimizations will include automatic frequency flattening.

from fireworks import Workflow
from atomate.qchem.fireworks.core import FrequencyFlatteningOptimizeFW, FrequencyFW
from atomate.utils.utils import get_logger

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "5/23/18"
__credits__ = "Brandon Wood, Shyam Dwaraknath"

logger = get_logger(__name__)


def get_wf_vac_pcm_smd(molecule,
                       pcm_dielectric,
                       smd_solvent,
                       linked=False,
                       qchem_input_params=None,
                       name="vac_pcm_smd",
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

    # Calculate frequencies with the vacuum optimized geometry and the pcm
    pcm_qchem_input_params = {"pcm_dielectric": pcm_dielectric}
    for key in orig_qchem_input_params:
        pcm_qchem_input_params[key] = orig_qchem_input_params[key]

    fw2 = FrequencyFW(
        name="{}:{}".format(molecule.composition.reduced_formula, "freq_pcm_vac_geom"+suffix),
        qchem_cmd=">>qchem_cmd<<",
        max_cores=">>max_cores<<",
        qchem_input_params=pcm_qchem_input_params,
        db_file=db_file,
        parents=fw1)

    # Calculate frequencies with the vacuum optimized geometry and the pcm
    smd_qchem_input_params = {"smd_solvent": smd_solvent}
    for key in orig_qchem_input_params:
        smd_qchem_input_params[key] = orig_qchem_input_params[key]

    fw3 = FrequencyFW(
        name="{}:{}".format(molecule.composition.reduced_formula, "freq_smd_vac_geom"+suffix),
        qchem_cmd=">>qchem_cmd<<",
        max_cores=">>max_cores<<",
        qchem_input_params=smd_qchem_input_params,
        db_file=db_file,
        parents=fw1)

    # Optimize in PCM
    fw4 = FrequencyFlatteningOptimizeFW(
        name="{}:{}".format(molecule.composition.reduced_formula, "pcm_FFopt"+suffix),
        qchem_cmd=">>qchem_cmd<<",
        max_cores=">>max_cores<<",
        qchem_input_params=pcm_qchem_input_params,
        linked=linked,
        db_file=db_file,
        parents=fw1)

    # Optimize in SMD
    fw5 = FrequencyFlatteningOptimizeFW(
        name="{}:{}".format(molecule.composition.reduced_formula, "smd_FFopt"+suffix),
        qchem_cmd=">>qchem_cmd<<",
        max_cores=">>max_cores<<",
        qchem_input_params=smd_qchem_input_params,
        linked=linked,
        db_file=db_file,
        parents=fw1)
      
    fws = [fw1, fw2, fw3, fw4, fw5]

    wfname = "{}:{}".format(molecule.composition.reduced_formula, name+suffix)

    return Workflow(fws, name=wfname, **kwargs)
