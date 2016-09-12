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

    Args:
        absorbing_atom:
        structure:
        radius:
        feff_input_set:
        feff_cmd:
        db_file:

    Returns:

    """
    fis = feff_input_set or MPXANESSet(absorbing_atom, structure, radius)
    fws = [XANESFW(absorbing_atom, structure, radius=radius, feff_input_set=fis, feff_cmd=feff_cmd,
                   db_file=db_file)]
    wfname = "{}:{}".format(structure.composition.reduced_formula, "XANES spectroscopy")
    return Workflow(fws, name=wfname)
