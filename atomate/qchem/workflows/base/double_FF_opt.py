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
__date__ = "5/23/18"
__credits__ = "Brandon Wood, Shyam Dwaraknath"

logger = get_logger(__name__)


def get_wf_double_FF_opt(molecule,
                         pcm_dielectric,
                         max_cores=32,
                         qchem_input_params=None,
                         name="douple_FF_opt",
                         qchem_cmd=">>qchem_cmd<<",
                         db_file=">>db_file<<",
                         **kwargs):
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
        max_cores (int): Maximum number of cores to parallelize over.
            Defaults to 32.
        qchem_input_params (dict): Specify kwargs for instantiating
            the input set parameters.
        qchem_cmd (str): Command to run QChem.
        db_file (str): path to file containing the database credentials.
        kwargs (keyword arguments): additional kwargs to be passed to Workflow

    Returns:
        Workflow
    """

    first_qchem_input_params = qchem_input_params or {}

    # Optimize the molecule in vacuum
    fw1 = FrequencyFlatteningOptimizeFW(
        molecule=molecule,
        name="first_FF_no_pcm",
        qchem_cmd=qchem_cmd,
        max_cores=max_cores,
        qchem_input_params=first_qchem_input_params,
        db_file=db_file)

    # Optimize the molecule in PCM
    second_qchem_input_params = {"pcm_dielectric": pcm_dielectric}
    for key in first_qchem_input_params:
        second_qchem_input_params[key] = first_qchem_input_params[key]
    fw2 = FrequencyFlatteningOptimizeFW(
        name="second_FF_with_pcm",
        qchem_cmd=qchem_cmd,
        max_cores=max_cores,
        qchem_input_params=second_qchem_input_params,
        db_file=db_file,
        parents=fw1)
    fws = [fw1, fw2]

    wfname = "{}:{}".format(molecule.composition.reduced_formula, name)

    return Workflow(fws, name=wfname, **kwargs)
