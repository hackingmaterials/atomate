# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines FEFF EXAFS spectroscopy workflows.
"""

from pymatgen.io.feff.sets import MPEXAFSSet

from fireworks import Workflow

from matmethods.utils.utils import get_logger
from matmethods.feff.fireworks.core import EXAFSFW

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


def get_wf_exafs(absorbing_atom, structure, radius=10.0, feff_input_set=None, feff_cmd="feff",
                 db_file=None, user_tag_settings=None):
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
    fis = feff_input_set or MPEXAFSSet(absorbing_atom, structure, radius,
                                       user_tag_settings=user_tag_settings or {})
    fws = [EXAFSFW(absorbing_atom, structure, radius=radius, feff_input_set=fis, feff_cmd=feff_cmd,
                   db_file=db_file)]
    print(fis.tags)
    wfname = "{}:{}".format(structure.composition.reduced_formula, "EXAFS spectroscopy")
    return Workflow(fws, name=wfname)
