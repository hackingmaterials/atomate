# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the gibbs free energy workflow.
"""

from datetime import datetime

from fireworks import Firework, Workflow

from pymatgen.analysis.elasticity.strain import Deformation

from matmethods.utils.utils import get_logger, append_fw_wf
from matmethods.vasp.firetasks.parse_outputs import GibbsFreeEnergyTask
from matmethods.vasp.workflows.base.deformations import get_wf_deformations

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


def get_wf_gibbs_free_energy(structure, vasp_input_set=None, vasp_cmd="vasp", deformations=None,
                             db_file=None, user_kpoints_settings=None, t_step=10, t_min=0, t_max=1000,
                             mesh=(20, 20, 20), eos="vinet"):
    """
    Returns quasi-harmonic gibbs free energy workflow.
    Note: phonopy package is required for the final analysis step.

    Args:
        structure (Structure): input structure.
        vasp_input_set (VaspInputSet)
        vasp_cmd (str): vasp command to run.
        deformations (list): list of deformation matrices(list of lists).
        db_file (str): path to the db file.
        user_kpoints_settings (dict): example: {"grid_density": 7000}
        t_step (float): temperature step (in K)
        t_min (float): min temperature (in K)
        t_max (float): max temperature (in K)
        mesh (list/tuple): reciprocal space density
        eos (str): equation of state used for fitting the energies and the volumes.
            supported equation of states: vinet, murnaghan, birch_murnaghan

    Returns:
        Workflow
    """
    try:
        from phonopy import Phonopy
    except ImportError:
        logger.warn("'phonopy' package NOT installed. It is required for the final analysis step.")

    tag = datetime.utcnow().strftime('%Y-%m-%d-%H-%M-%S-%f')

    deformations = [Deformation(defo_mat) for defo_mat in deformations]
    wf_gibbs = get_wf_deformations(structure, deformations, name="gibbs deformation",
                                   vasp_input_set=vasp_input_set, lepsilon=True, vasp_cmd=vasp_cmd,
                                   db_file=db_file, user_kpoints_settings=user_kpoints_settings,
                                   tag=tag)

    fw_analysis = Firework(GibbsFreeEnergyTask(tag=tag, db_file=db_file, t_step=t_step, t_min=t_min,
                                               t_max=t_max, mesh=mesh, eos=eos),
                           name="Gibbs Free Energy")

    append_fw_wf(wf_gibbs, fw_analysis)

    wf_gibbs.name = "{}:{}".format(structure.composition.reduced_formula, "gibbs free energy")

    return wf_gibbs
