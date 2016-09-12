# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines FEFF XANES spectroscopy workflows.
"""

from pymatgen.io.feff.sets import MPXANESSet

from fireworks import Workflow

from matmethods.utils.utils import get_logger
from matmethods.feff.fireworks.core import XANESFW

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


def get_wf_xanes(absorbing_atom, structure, radius=10.0, feff_input_set=None, feff_cmd="feff",
                 db_file=None):
    """
    Returns FEFF XANES spectroscopy workflow.

    Args:
        absorbing_atom (str): absorbing atom symbol
        structure (Structure): input structure
        radius (float): cluster radius in angstroms
        feff_input_set (FeffDictSet): the input set for the FEFF run
        feff_cmd (str): path to the feff binary
        db_file (str):  path to the db file.

    Returns:
        Workflow
    """
    fis = feff_input_set or MPXANESSet(absorbing_atom, structure, radius)
    fws = [XANESFW(absorbing_atom, structure, radius=radius, feff_input_set=fis, feff_cmd=feff_cmd,
                   db_file=db_file)]
    wfname = "{}:{}".format(structure.composition.reduced_formula, "XANES spectroscopy")
    return Workflow(fws, name=wfname)
