# coding: utf-8

# This module defines a workflow that fragments a molecule and optimizes each fragment.
# It will most often be used in order to obtain the bond dissociation energies of a molecule.

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


def get_fragmentation_wf(molecule,
                         depth=1,
                         open_rings=True,
                         additional_charges=None,
                         do_triplets=True,
                         pcm_dielectric=None,
                         do_optimization=True,
                         linked=False,
                         qchem_input_params=None,
                         name="FF then fragment",
                         db_file=">>db_file<<",
                         check_db=True,
                         **kwargs):
    """

    Args:
        molecule (Molecule): input molecule to be fragmented.
        depth (int): The number of levels of iterative fragmentation to perform,
                     where each evel will include fragments obtained by breaking
                     one bond of a fragment one level up. If set to 0, instead
                     all possible fragments are generated using an alternative,
                     non-iterative scheme. Defaults to 1.
        open_rings (bool): Whether or not to open any rings encountered during fragmentation.
                           Defaults to True. If true, any bond that fails to yield disconnected
                           graphs when broken is instead removed and the entire structure is
                           optimized with OpenBabel in order to obtain a good initial guess for
                           an opened geometry that can then be put back into QChem to be
                           optimized without the ring just reforming.
        additional_charges (list): List of additional charges besides the defaults described in the
                                   firetask. For example, if a principle molecule with a +2 charge
                                   is provided, by default all fragments will be calculated with
                                   +1 and +2 charges. If the user includes additional_charges=[0]
                                   then all fragments will be calculated with 0, +1, and +2 charges.
                                   Additional charge values of 1 or 2 would not cause any new charges
                                   to be calculated as they are already done. Defaults to [].
        do_triplets (bool): Whether to simulate triplets as well as singlets for molecules with
                            an even number of electrons. Defaults to True.
        pcm_dielectric (float): The PCM dielectric constant.
        do_optimization (bool): Whether or not to optimize the given molecule
                                before fragmentation. Defaults to True.
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
        check_db (bool): Whether or not to check the database for equivalent
                         structures before adding new fragment fireworks.
                         Defaults to True.
        kwargs (keyword arguments): additional kwargs to be passed to Workflow

    Returns:
        Workflow with the following fireworks:

        Firework 1 : write QChem input for an FF optimization,
                     run FF_opt QCJob,
                     parse directory and insert into db,
                     pass relaxed molecule to fw_spec and on to fw2,

        Firework 2 : find all unique fragments of the optimized molecule
                     and add a frequency flattening optimize FW to the
                     workflow for each one

        Note that Firework 1 is only present if do_optimization=True.

    """

    qchem_input_params = qchem_input_params or {}
    additional_charges = additional_charges or []
    if pcm_dielectric != None:
        qchem_input_params["pcm_dielectric"] = pcm_dielectric

    if do_optimization:
        # Optimize the original molecule
        fw1 = FrequencyFlatteningOptimizeFW(
            molecule=molecule,
            name="first FF",
            qchem_cmd=">>qchem_cmd<<",
            max_cores=">>max_cores<<",
            qchem_input_params=qchem_input_params,
            linked=linked,
            db_file=db_file)

        # Fragment the optimized molecule
        fw2 = FragmentFW(
            depth=depth,
            open_rings=open_rings,
            additional_charges=additional_charges,
            do_triplets=do_triplets,
            linked=linked,
            name="fragment and FF_opt",
            qchem_input_params=qchem_input_params,
            db_file=db_file,
            check_db=check_db,
            parents=fw1)
        fws = [fw1, fw2]

    else:
        # Fragment the given molecule
        fw1 = FragmentFW(
            molecule=molecule,
            depth=depth,
            open_rings=open_rings,
            additional_charges=additional_charges,
            do_triplets=do_triplets,
            linked=linked,
            name="fragment and FF_opt",
            qchem_input_params=qchem_input_params,
            db_file=db_file,
            check_db=check_db)
        fws = [fw1]

    wfname = "{}:{}".format(molecule.composition.reduced_formula, name)

    return Workflow(fws, name=wfname, **kwargs)
