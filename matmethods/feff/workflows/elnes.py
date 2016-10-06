# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines FEFF ELNES spectroscopy workflows.
"""

from pymatgen.util.testing import PymatgenTest
from pymatgen.io.feff.sets import MPELNESSet

from fireworks import Workflow

from matmethods.utils.utils import get_logger
from matmethods.feff.fireworks.core import ELNESFW

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


def get_wf_elnes(absorbing_atom, structure=None, edge="K", radius=10., beam_energy=100,
                 beam_direction=None, collection_angle=1, convergence_angle=1,
                 user_eels_settings=None, feff_input_set=None, feff_cmd="feff", db_file=None,
                 mp_id=None):
    """
    Returns FEFF-ELNES spectroscopy workflow.

    Args:
        absorbing_atom (str): absorbing atom symbol
        structure (Structure): input structure. If None and mp_id is provided, the corresponding
            structure will be fetched from the Materials Project db.
        edge (str): absorption edge
        radius (float): cluster radius in angstroms
        beam_energy (float): the incident beam energy in keV
        beam_direction (list): incident beam direction. Default is none ==> the spectrum will be
            averaged over all directions.
        collection_angle (float): collection angle in mrad
        convergence_angle (float): convergence angle in mrad
        user_eels_settings (dict): override default eels settings, has effect only when
            feff_input_set = None
        feff_input_set (FeffDictSet): the input set for the FEFF run
        feff_cmd (str): path to the feff binary
        db_file (str):  path to the db file.
        mp_id (str): mp id of the input structure.

    Returns:
        Workflow
    """
    if structure:
        wfname = "{}:{}".format(structure.composition.reduced_formula, "ELNES spectroscopy")
    elif mp_id:
            structure = PymatgenTest.get_mp_structure(mp_id)
            wfname = "{}-{}:{}".format(mp_id, structure.composition.reduced_formula,
                                       "ELNES spectroscopy")
    else:
        raise ValueError("Neither structure no mp_id provided.")
    fis = feff_input_set or MPELNESSet(absorbing_atom, structure, edge, radius, beam_energy,
                                       beam_direction, collection_angle, convergence_angle,
                                       user_eels_settings=user_eels_settings)
    fws = [ELNESFW(absorbing_atom, structure, edge=edge, radius=radius, beam_energy=beam_energy,
                   beam_direction=beam_direction, collection_angle=collection_angle,
                   convergence_angle=convergence_angle, user_eels_settings=user_eels_settings,
                   feff_input_set=fis, feff_cmd=feff_cmd, db_file=db_file, mp_id=mp_id)]
    return Workflow(fws, name=wfname)
