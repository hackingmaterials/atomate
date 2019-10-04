# coding: utf-8


"""
This module defines the gibbs free energy workflow.
"""

from uuid import uuid4

from fireworks import Firework, Workflow

from pymatgen.analysis.elasticity.strain import Deformation
from pymatgen.io.vasp.sets import MPStaticSet

from atomate.utils.utils import get_logger
from atomate.vasp.firetasks.parse_outputs import GibbsAnalysisToDb
from atomate.vasp.workflows.base.deformations import get_wf_deformations

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


def get_wf_gibbs_free_energy(structure, deformations, vasp_input_set=None, vasp_cmd="vasp",
                             db_file=None, user_kpoints_settings=None, t_step=10, t_min=0,
                             t_max=1000, mesh=(20, 20, 20), eos="vinet", qha_type="debye_model",
                             pressure=0.0, poisson=0.25, anharmonic_contribution=False,
                             metadata=None, tag=None):
    """
    Returns quasi-harmonic gibbs free energy workflow.
    Note: phonopy package is required for the final analysis step if qha_type="phonopy"

    Args:
        structure (Structure): input structure.
        deformations (list): list of deformation matrices(list of lists).
        vasp_input_set (VaspInputSet)
        vasp_cmd (str): vasp command to run.
        db_file (str): path to the db file.
        user_kpoints_settings (dict): example: {"grid_density": 7000}
        t_step (float): temperature step (in K)
        t_min (float): min temperature (in K)
        t_max (float): max temperature (in K)
        mesh (list/tuple): reciprocal space density
        eos (str): equation of state used for fitting the energies and the volumes.
            options supported by phonopy: "vinet", "murnaghan", "birch_murnaghan".
            Note: pymatgen supports more options than phonopy. see pymatgen.analysis.eos.py
        qha_type(str): quasi-harmonic approximation type: "debye_model" or "phonopy",
            default is "debye_model"
        pressure (float): in GPa
        poisson (float): poisson ratio
        anharmonic_contribution (bool): consider anharmonic contributions to
            Gibbs energy from the Debye model. Defaults to False.
        metadata (dict): meta data
        tag (str): something unique to identify the tasks in this workflow. If None a random uuid
            will be assigned.

    Returns:
        Workflow
    """

    tag = tag or "gibbs group: >>{}<<".format(str(uuid4()))

    deformations = [Deformation(defo_mat) for defo_mat in deformations]

    # static input set for the transmuter fireworks
    vis_static = vasp_input_set
    if vis_static is None:
        lepsilon = False
        if qha_type not in ["debye_model"]:
            lepsilon = True
            try:
                from phonopy import Phonopy
            except ImportError:
                raise RuntimeError("'phonopy' package is NOT installed but is required for the final "
                                   "analysis step; you can alternatively switch to the qha_type to "
                                   "'debye_model' which does not require 'phonopy'.")
        vis_static = MPStaticSet(structure, force_gamma=True, lepsilon=lepsilon,
                                 user_kpoints_settings=user_kpoints_settings)

    wf_gibbs = get_wf_deformations(structure, deformations, name="gibbs deformation",
                                   vasp_cmd=vasp_cmd, db_file=db_file, tag=tag, metadata=metadata,
                                   vasp_input_set=vis_static)

    fw_analysis = Firework(GibbsAnalysisToDb(tag=tag, db_file=db_file, t_step=t_step, t_min=t_min,
                                             t_max=t_max, mesh=mesh, eos=eos, qha_type=qha_type,
                                             pressure=pressure, poisson=poisson, metadata=metadata,
                                             anharmonic_contribution=anharmonic_contribution,),
                           name="Gibbs Free Energy")

    wf_gibbs.append_wf(Workflow.from_Firework(fw_analysis), wf_gibbs.leaf_fw_ids)

    wf_gibbs.name = "{}:{}".format(structure.composition.reduced_formula, "gibbs free energy")

    return wf_gibbs
