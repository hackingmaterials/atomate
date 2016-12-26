# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the bulk modulus workflow.
"""

from datetime import datetime

from fireworks import Firework, Workflow

from pymatgen.analysis.elasticity.strain import Deformation

from atomate.utils.utils import get_logger, append_fw_wf
from atomate.vasp.firetasks.parse_outputs import FitEquationOfStateTask
from atomate.vasp.workflows.base.deformations import get_wf_deformations

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


def get_wf_bulk_modulus(structure, deformations, vasp_input_set=None, vasp_cmd="vasp", db_file=None,
                        user_kpoints_settings=None, eos="vinet"):
    """
    Returns the workflow that computes the bulk modulus by fitting to the given equation of state.

    Args:
        structure (Structure): input structure.
        deformations (list): list of deformation matrices(list of lists).
        vasp_input_set (VaspInputSet)
        vasp_cmd (str): vasp command to run.
        db_file (str): path to the db file.
        user_kpoints_settings (dict): example: {"grid_density": 7000}
        eos (str): equation of state used for fitting the energies and the volumes.
            supported equation of states: "quadratic", "murnaghan", "birch", "birch_murnaghan",
            "pourier_tarantola", "vinet", "deltafactor". See pymatgen.analysis.eos.py

    Returns:
        Workflow
    """

    tag = datetime.utcnow().strftime('%Y-%m-%d-%H-%M-%S-%f')

    deformations = [Deformation(defo_mat) for defo_mat in deformations]
    wf_bulk_modulus = get_wf_deformations(structure, deformations, name="bulk_modulus deformation",
                                          vasp_input_set=vasp_input_set, lepsilon=False,
                                          vasp_cmd=vasp_cmd, db_file=db_file,
                                          user_kpoints_settings=user_kpoints_settings, tag=tag)

    fw_analysis = Firework(FitEquationOfStateTask(tag=tag, db_file=db_file, eos=eos),
                           name="fit equation of state")

    append_fw_wf(wf_bulk_modulus, fw_analysis)

    wf_bulk_modulus.name = "{}:{}".format(structure.composition.reduced_formula, "Bulk modulus")

    return wf_bulk_modulus
