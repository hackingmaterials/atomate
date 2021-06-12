# This module defines a workflow for FFopting a molecule and then analyzing its
# electron density critical points with Critic2.

from fireworks import Workflow
from atomate.qchem.fireworks.core import FrequencyFlatteningOptimizeFW, CubeAndCritic2FW
from atomate.utils.utils import get_logger

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "11/20/19"

logger = get_logger(__name__)


def get_wf_FFopt_and_critic(
    molecule, suffix, qchem_input_params=None, db_file=">>db_file<<", **kwargs
):
    """

    Firework 1 : write QChem input for an FF optimization,
                 run FF_opt QCJob,
                 parse directory and insert into db,
                 pass relaxed molecule to fw_spec and on to fw2,

    Firework 2 : write QChem input for a single point calc to print a cube file
                 run SP QCJob, thereby printing a cube file
                 run Critic2 on the printed cube file
                 parse directory and insert into db

    Args:
        molecule (Molecule): input molecule to be optimized and run.
        suffix (str): Workflow naming suffix
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
        db_file (str): path to file containing the database credentials.
        kwargs (keyword arguments): additional kwargs to be passed to Workflow

    Returns:
        Workflow
    """

    qchem_input_params = qchem_input_params or {}

    # FFopt
    fw1 = FrequencyFlatteningOptimizeFW(
        molecule=molecule,
        name="{}:{}".format(
            molecule.composition.alphabetical_formula, "FFopt_" + suffix
        ),
        qchem_cmd=">>qchem_cmd<<",
        max_cores=">>max_cores<<",
        qchem_input_params=qchem_input_params,
        linked=True,
        db_file=">>db_file<<",
    )

    # Critic
    fw2 = CubeAndCritic2FW(
        name="{}:{}".format(molecule.composition.alphabetical_formula, "CC2_" + suffix),
        qchem_cmd=">>qchem_cmd<<",
        max_cores=">>max_cores<<",
        qchem_input_params=qchem_input_params,
        db_file=">>db_file<<",
        parents=fw1,
    )
    fws = [fw1, fw2]

    wfname = "{}:{}".format(
        molecule.composition.alphabetical_formula, "FFopt_CC2_WF_" + suffix
    )

    return Workflow(fws, name=wfname, **kwargs)
