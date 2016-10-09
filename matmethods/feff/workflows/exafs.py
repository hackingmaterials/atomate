# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines FEFF EXAFS spectroscopy workflow.
"""

from pymatgen.io.feff.sets import MPEXAFSSet

from fireworks import Workflow

from matmethods.utils.utils import get_logger
from matmethods.feff.fireworks.core import EXAFSFW

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


def get_wf_exafs(absorbing_atom, structure, edge="K", radius=10.0, feff_input_set=None,
                 feff_cmd="feff", db_file=None, user_tag_settings=None):
    """
    Returns FEFF EXAFS spectroscopy workflow.


    Args:
        absorbing_atom (str): absorbing atom symbol
        structure (Structure): input structure
        edge (str): absorption edge
        radius (float): cluster radius in angstroms
        feff_input_set (FeffDictSet): the input set for the FEFF run
        feff_cmd (str): path to the feff binary
        db_file (str):  path to the db file.

    Returns:
        Workflow
    """
    fis = feff_input_set or MPEXAFSSet(absorbing_atom, structure, edge=edge, radius=radius,
                                       user_tag_settings=user_tag_settings or {})
    fws = [EXAFSFW(absorbing_atom, structure, edge=edge, radius=radius, feff_input_set=fis,
                   feff_cmd=feff_cmd, db_file=db_file)]

    wfname = "{}:{}:{} edge".format(structure.composition.reduced_formula, "EXAFS spectroscopy", edge)
    return Workflow(fws, name=wfname)
