# coding: utf-8

from __future__ import absolute_import, division, print_function, \
    unicode_literals

# This module defines a workflow that first performs a frequency flattening
# optimization of a molecule in vacuum and then finds all unique fragments
# of that molecule and performs a frequency flattening optimization on each.

from fireworks import Workflow
from atomate.qchem.fireworks.core import FrequencyFlatteningOptimizeFW, FragmentFW
from atomate.utils.utils import get_logger

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "6/22/18"
__credits__ = "Brandon Wood, Shyam Dwaraknath"

logger = get_logger(__name__)


def get_wf_FF_then_fragment(molecule,
                            depth,
                            pcm_dielectric=None,
                            max_cores=">>max_cores<<",
                            qchem_input_params=None,
                            name="FF then fragment",
                            qchem_cmd=">>qchem_cmd<<",
                            db_file=">>db_file<<",
                            check_db=True,
                            **kwargs):
    """

    Firework 1 : write QChem input for an FF optimization,
                 run FF_opt QCJob,
                 parse directory and insert into db,
                 pass relaxed molecule to fw_spec and on to fw2,

    Firework 2 : find all unique fragments of the optimized molecule
                 and add a frequency flattening optimize FW to the
                 workflow for each one

    Args:
        molecule (Molecule): input molecule to be optimized and run.
        pcm_dielectric (float): The PCM dielectric constant.
        max_cores (int): Maximum number of cores to parallelize over.
            Defaults to 32.
        qchem_input_params (dict): Specify kwargs for instantiating
            the input set parameters.
        qchem_cmd (str): Command to run QChem.
        db_file (str): path to file containing the database credentials.
        check_db (bool): Whether or not to check the database for equivalent
            structures before adding new fragment fireworks. Defaults to True.
        kwargs (keyword arguments): additional kwargs to be passed to Workflow

    Returns:
        Workflow
    """

    qchem_input_params = qchem_input_params or {}
    if pcm_dielectric != None:
        qchem_input_params["pcm_dielectric"] = pcm_dielectric

    # Optimize the original molecule
    fw1 = FrequencyFlatteningOptimizeFW(
        molecule=molecule,
        name="first FF",
        qchem_cmd=qchem_cmd,
        max_cores=max_cores,
        qchem_input_params=qchem_input_params,
        db_file=db_file)

    # Fragment the optimized molecule
    fw2 = FragmentFW(
        depth=depth,
        name="fragment and FF_opt",
        qchem_cmd=qchem_cmd,
        max_cores=max_cores,
        qchem_input_params=qchem_input_params,
        db_file=db_file,
        check_db=check_db,
        parents=fw1)
    fws = [fw1, fw2]

    wfname = "{}:{}".format(molecule.composition.reduced_formula, name)

    return Workflow(fws, name=wfname, **kwargs)
