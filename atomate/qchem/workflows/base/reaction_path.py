# coding: utf-8

from __future__ import absolute_import, division, print_function, \
    unicode_literals

# This module defines a workflow for optimizing a transition-state geometry
# and then identifying the reaction path

from fireworks import Workflow
from atomate.qchem.fireworks.core import (FrequencyFlatteningOptimizeFW)
from atomate.utils.utils import get_logger

__author__ = "Evan Spotte-Smith"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Alpha"
__date__ = "February 2020"
__credits__ = "Sam Blau, Alex Epstein, Trevor Seguin"

logger = get_logger(__name__)


def get_wf_reaction_path_with_ts(molecule,
                                 mode,
                                 suffix,
                                 scale=1.0,
                                 qchem_cmd=">>qchem_cmd<<",
                                 max_cores=">>max_cores<<",
                                 multimode=">>multimode<<",
                                 qchem_input_params=None,
                                 name="reaction_path_with_ts",
                                 db_file=">>db_file<<",
                                 **kwargs):
    """
    Firework 1: Perturb the given molecule along the given frequency
                mode in the forwards direction,
                run FF_opt QCJob,
                parse directory, and insert into db
    Firework 2: Perturb the given molecule along the given frequency
                mode in the reverse direction,
                run FF_opt QCJob,
                parse directory, and insert into db
    Args:
        molecule (Molecule): pymatgen Molecule representing the TS.
        mode (np.ndarray): a matrix used to perturb the molecule
        scale (float): Scaling factor for molecular perturbations.
        qchem_cmd (str): Command to run QChem.
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
        name (string): name for the Workflow
        db_file (str): path to file containing the database credentials.
        kwargs (keyword arguments): additional kwargs to be passed to Workflow
    Returns:
        Workflow
    """

    fw1 = FrequencyFlatteningOptimizeFW(
        molecule=molecule,
        name="{}:perturb_forwards{}".format(molecule.composition.alphabetical_formula, suffix),
        qchem_cmd=qchem_cmd,
        max_cores=max_cores,
        multimode=multimode,
        qchem_input_params=qchem_input_params,
        perturb_geometry=True,
        scale=1.0 * scale,
        mode=mode,
        linked=True,
        db_file=db_file)

    fw2 = FrequencyFlatteningOptimizeFW(
        molecule=molecule,
        name="{}:perturb_backwards{}".format(molecule.composition.alphabetical_formula, suffix),
        qchem_cmd=qchem_cmd,
        max_cores=max_cores,
        multimode=multimode,
        qchem_input_params=qchem_input_params,
        perturb_geometry=True,
        scale=-1.0 * scale,
        mode=mode,
        linked=True,
        db_file=db_file)

    fws = [fw1, fw2]

    wfname = "{}:{}".format(molecule.composition.alphabetical_formula, name)

    return Workflow(fws, name=wfname, **kwargs)