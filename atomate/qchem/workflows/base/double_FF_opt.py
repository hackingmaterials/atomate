# coding: utf-8

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
                         linked=False,
                         qchem_input_params=None,
                         name="douple_FF_opt",
                         db_file=">>db_file<<",
                         **kwargs):
    """

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
        qchem_input_params (dict): Specify kwargs for instantiating the input set parameters.
                                   Basic uses would be to modify the default inputs of the set,
                                   such as dft_rung, basis_set, pcm_dielectric, scf_algorithm,
                                   or max_scf_cycles. See pymatgen/io/qchem/sets.py for default
                                   values of all input parameters. For instance, if a user wanted
                                   to use a more advanced DFT functional, include a pcm with a
                                   dielectric of 30, and use a larger basis, the user would set
                                   qchem_input_params = {"dft_rung": 5, "pcm_dielectric": 30,
                                   "basis_set": "6-311++g**"}. However, more advanced customization
                                   of the input is also possible through the overwrite_inputs key
                                   which allows the user to directly modify the rem, pcm, smd, and
                                   solvent dictionaries that QChemDictSet passes to inputs.py to
                                   print an actual input file. For instance, if a user wanted to
                                   set the sym_ignore flag in the rem section of the input file
                                   to true, then they would set qchem_input_params = {"overwrite_inputs":
                                   "rem": {"sym_ignore": "true"}}. Of course, overwrite_inputs
                                   could be used in conjuction with more typical modifications,
                                   as seen in the test_double_FF_opt workflow test.
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
        qchem_cmd=">>qchem_cmd<<",
        max_cores=">>max_cores<<",
        qchem_input_params=first_qchem_input_params,
        linked=linked,
        db_file=db_file)

    # Optimize the molecule in PCM
    second_qchem_input_params = {"pcm_dielectric": pcm_dielectric}
    for key in first_qchem_input_params:
        second_qchem_input_params[key] = first_qchem_input_params[key]
    fw2 = FrequencyFlatteningOptimizeFW(
        name="second_FF_with_pcm",
        qchem_cmd=">>qchem_cmd<<",
        max_cores=">>max_cores<<",
        qchem_input_params=second_qchem_input_params,
        linked=linked,
        db_file=db_file,
        parents=fw1)
    fws = [fw1, fw2]

    wfname = "{}:{}".format(molecule.composition.reduced_formula, name)

    return Workflow(fws, name=wfname, **kwargs)
