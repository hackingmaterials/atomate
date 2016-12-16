# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the thermal expansion coefficient workflow.
"""

from datetime import datetime

from fireworks import Firework, Workflow

from pymatgen.analysis.elasticity.strain import Deformation

from atomate.utils.utils import get_logger, append_fw_wf
from atomate.vasp.firetasks.parse_outputs import ThermalExpansionCoeffTask
from atomate.vasp.workflows.base.deformations import get_wf_deformations

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


def get_wf_thermal_expansion(structure, deformations, vasp_input_set=None, vasp_cmd="vasp",
                             db_file=None, user_kpoints_settings=None, t_step=10, t_min=0, t_max=1000,
                             mesh=(20, 20, 20), eos="vinet", pressure=0.0):
    """
    Returns quasi-harmonic thermal expansion workflow.
    Note: phonopy package is required for the final analysis step.

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
        pressure (float): in GPa

    Returns:
        Workflow
    """
    try:
        from phonopy import Phonopy
    except ImportError:
        logger.warn("'phonopy' package NOT installed. Required for the final analysis step.")

    tag = datetime.utcnow().strftime('%Y-%m-%d-%H-%M-%S-%f')

    deformations = [Deformation(defo_mat) for defo_mat in deformations]
    wf_alpha = get_wf_deformations(structure, deformations, name="thermal_expansion deformation",
                                   vasp_input_set=vasp_input_set, lepsilon=True, vasp_cmd=vasp_cmd,
                                   db_file=db_file, user_kpoints_settings=user_kpoints_settings,
                                   tag=tag)

    fw_analysis = Firework(ThermalExpansionCoeffTask(tag=tag, db_file=db_file, t_step=t_step,
                                                     t_min=t_min, t_max=t_max, mesh=mesh, eos=eos,
                                                     pressure=pressure),
                           name="Thermal expansion")

    append_fw_wf(wf_alpha, fw_analysis)

    wf_alpha.name = "{}:{}".format(structure.composition.reduced_formula, "thermal expansion")

    return wf_alpha
