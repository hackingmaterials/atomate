# coding: utf-8

from __future__ import absolute_import, division, print_function, \
    unicode_literals

"""
This module defines a workflow for optimizing a molecule first in vacuum and then
in PCM. Both optimizations will include automatic frequency flattening. 
"""

import numpy as np

from fireworks import Firework, Workflow

from atomate.qchem.fireworks.core import FrequencyFlatteningOptimizeFW
from atomate.utils.utils import get_logger

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"

logger = get_logger(__name__)


def get_wf_double_FF_opt(molecule, pcm_dielectric, name="douple_FF_opt", db_file=None, **kwargs):
    """
    Returns a workflow to the torsion potential for a molecule.

    Firework 1 : write QChem input for an FF optimization,
                 run FF_opt QCJob,
                 parse directory and insert into db,
                 pass relaxed molecule to fw_spec and on to fw2,

    Firework 2 : write QChem input for an optimization in the 
                    presence of a PCM, using the molecule passed
                    from fw1,
                 run FF_opt QCJob,
                 parse directory and insert into db

    Args:
        molecule (Molecule): input molecule to be optimized and run.
        pcm_dielectric (float): The PCM dielectric constant.
        db_file (str): path to file containing the database credentials.
        kwargs (keyword arguments): additional kwargs to be passed to Workflow

    Returns:
        Workflow
    """

    # Optimize the molecule in vacuum
    fw1 = FrequencyFlatteningOptimizeFW(molecule=molecule)
    # Optimize the molecule in PCM
    fw2 = FrequencyFlatteningOptimizeFW(parents=fw1, qchem_input_params={"pcm_dielectric": pcm_dielectric})
    fws = [fw1, fw2]

    wfname = "{}:{}".format(Molecule.formula_name, name)

    return Workflow(fws, name=wfname, **kwargs)
    


