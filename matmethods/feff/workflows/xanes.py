# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines FEFF XANES spectroscopy workflows.
"""

from pymatgen.io.feff.sets import MPXANESSet

from fireworks import Workflow

from matmethods.utils.utils import get_logger
from matmethods.feff.fireworks.core import XASFW
from matmethods.feff.utils import get_all_absorbing_atoms


__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


def get_wf_xanes(absorbing_atom, structure, edge="K", radius=10.0, feff_input_set=None,
                 feff_cmd="feff", db_file=None, metadata=None):
    """
    Returns FEFF XANES spectroscopy workflow.

    Args:
        absorbing_atom (str/int): absorbing atom symbol or site index
        structure (Structure): input structure
        edge (str): absorption edge
        radius (float): cluster radius in angstroms
        feff_input_set (FeffDictSet): the input set for the FEFF run
        feff_cmd (str): path to the feff binary
        db_file (str):  path to the db file.
        metadata (dict): meta data

    Returns:
        Workflow
    """
    structure = structure.get_primitive_structure()

    # get the absorbing atom site index/indices
    ab_atom_indices = get_all_absorbing_atoms(absorbing_atom, structure)

    fws = []
    for ab_idx in ab_atom_indices:
        fw_metadata = dict(metadata) if metadata else {}
        fw_metadata["absorbing_atom_index"] = ab_idx
        fis = feff_input_set or MPXANESSet(absorbing_atom, structure, edge=edge, radius=radius)
        fws.append(XASFW(absorbing_atom, structure, "XANES", edge=edge, radius=radius,
                         feff_input_set=fis, feff_cmd=feff_cmd, db_file=db_file, metadata=fw_metadata,
                         name="XANES spectroscopy"))

    wf_metadata = dict(metadata) if metadata else {}
    wf_metadata["absorbing_atom_indices"] = list(ab_atom_indices)
    wfname = "{}:{}:{} edge".format(structure.composition.reduced_formula, "XANES spectroscopy", edge)

    return Workflow(fws, name=wfname, metadata=metadata)
