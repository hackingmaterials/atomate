# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the gibbs free energy workflow.
"""

from uuid import uuid4

from fireworks import Firework, Workflow

from pymatgen.analysis.elasticity.strain import Deformation

from atomate.utils.utils import get_logger
from atomate.vasp.firetasks.parse_outputs import GibbsAnalysisToDb
from atomate.vasp.workflows.base.deformations import get_wf_deformations

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


def get_wf_gibbs_free_energy(structure, deformations, vasp_input_set=None, vasp_cmd="vasp",
                             db_file=None, user_kpoints_settings=None, t_step=10, t_min=0, t_max=1000,
                             mesh=(20, 20, 20), eos="vinet", qha_type="debye_model", pressure=0.0,
                             poisson=0.25, metadata=None):
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
        metadata (dict): meta data

    Returns:
        Workflow
    """
    lepsilon = False
    if qha_type not in ["debye_model"]:
        lepsilon = True
        try:
            from phonopy import Phonopy
        except ImportError:
            raise RuntimeError("'phonopy' package is NOT installed but is required for the final "
                               "analysis step; you can alternatively switch to the qha_type to "
                               "'debye_model' which does not require 'phonopy'.")

    tag = "gibbs group: >>{}<<".format(str(uuid4()))

    deformations = [Deformation(defo_mat) for defo_mat in deformations]
    wf_gibbs = get_wf_deformations(structure, deformations, name="gibbs deformation",
                                   vasp_input_set=vasp_input_set, lepsilon=lepsilon, vasp_cmd=vasp_cmd,
                                   db_file=db_file, user_kpoints_settings=user_kpoints_settings,
                                   tag=tag, metadata=metadata)

    # TODO: @kmathew - better to stick to lowercase FW names ('gibbs free energy'). That is the
    # convention for most FWs and switching back and forth is confusing to remember. -computron
    fw_analysis = Firework(GibbsAnalysisToDb(tag=tag, db_file=db_file, t_step=t_step, t_min=t_min,
                                             t_max=t_max, mesh=mesh, eos=eos, qha_type=qha_type,
                                             pressure=pressure, poisson=poisson, metadata=metadata),
                           name="Gibbs Free Energy")

    wf_gibbs.append_wf(Workflow.from_Firework(fw_analysis), wf_gibbs.leaf_fw_ids)

    wf_gibbs.name = "{}:{}".format(structure.composition.reduced_formula, "gibbs free energy")

    return wf_gibbs
