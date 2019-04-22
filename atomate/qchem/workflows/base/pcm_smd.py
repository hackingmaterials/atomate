# coding: utf-8

from __future__ import absolute_import, division, print_function, \
    unicode_literals

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


def get_wf_pcm_smd(molecule,
                   pcm_dielectric,
                   smd_solvent,
                   linked=True,
                   qchem_input_params=None,
                   name="pcm_smd",
                   suffix="",
                   db_file=">>db_file<<",
                   **kwargs):
    """

    Returns:
        Workflow
    """

    orig_qchem_input_params = qchem_input_params or {}
    if smd_solvent == "custom" or smd_solvent == "other":
        if "custom_smd" not in orig_qchem_input_params:
            raise RuntimeError(
                'Custom SMD parameters must be provided in qchem_input_params' +
                '["custom_smd"] as a string of seven comma separated values in the ' +
                'following order: dielectric, refractive index, acidity, basicity, ' +
                'surface tension, aromaticity, electronegative halogenicity'
                )

    # Calculate frequencies with the vacuum optimized geometry and the pcm
    pcm_qchem_input_params = {"pcm_dielectric": pcm_dielectric}
    for key in orig_qchem_input_params:
        pcm_qchem_input_params[key] = orig_qchem_input_params[key]

    # Calculate frequencies with the vacuum optimized geometry and the pcm
    smd_qchem_input_params = {"smd_solvent": smd_solvent}
    for key in orig_qchem_input_params:
        smd_qchem_input_params[key] = orig_qchem_input_params[key]

    # Optimize in PCM
    fw1 = FrequencyFlatteningOptimizeFW(
        molecule=molecule,
        name="{}:{}".format(molecule.composition.reduced_formula, "pcm_FFopt"+suffix),
        qchem_cmd=">>qchem_cmd<<",
        max_cores=">>max_cores<<",
        qchem_input_params=pcm_qchem_input_params,
        linked=linked,
        db_file=db_file)

    # Optimize in SMD
    fw2 = FrequencyFlatteningOptimizeFW(
        molecule=molecule,
        name="{}:{}".format(molecule.composition.reduced_formula, "smd_FFopt"+suffix),
        qchem_cmd=">>qchem_cmd<<",
        max_cores=">>max_cores<<",
        qchem_input_params=smd_qchem_input_params,
        linked=linked,
        db_file=db_file)
      
    fws = [fw1, fw2]

    wfname = "{}:{}".format(molecule.composition.reduced_formula, name+suffix)

    return Workflow(fws, name=wfname, **kwargs)
