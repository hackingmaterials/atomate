# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines FEFF EELS(ELNES/EXELFS) spectroscopy workflows.
"""

from fireworks import Workflow

from atomate.utils.utils import get_logger
from atomate.feff.fireworks.core import EELSFW
from atomate.feff.utils import get_all_absorbing_atoms

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


def get_wf_eels(absorbing_atom, structure=None, spectrum_type="ELNES", edge="K", radius=10.,
                beam_energy=100, beam_direction=None, collection_angle=1, convergence_angle=1,
                user_eels_settings=None, user_tag_settings=None, feff_cmd="feff", db_file=None,
                metadata=None, use_primitive=False, feff_input_set=None):
    """
    Returns FEFF ELNES/EXELFS spectroscopy workflow.

    Args:
        absorbing_atom (str): absorbing atom symbol
        structure (Structure): input structure. If None and mp_id is provided, the corresponding
            structure will be fetched from the Materials Project db.
        spectrum_type (str): ELNES or EXELFS
        edge (str): absorption edge. K, L1, L2, L3
        radius (float): cluster radius in angstroms. Ignored for reciprocal space calculations
        beam_energy (float): the incident beam energy in keV
        beam_direction (list): incident beam direction. Default is none ==> the spectrum will be
            averaged over all directions.
        collection_angle (float): collection angle in mrad
        convergence_angle (float): convergence angle in mrad
        user_eels_settings (dict): override default eels settings.
        user_tag_settings (dict): override other general feff default tag settings.
        feff_cmd (str): path to the feff binary
        db_file (str):  path to the db file.
        metadata (dict): meta data
        use_primitive (bool): convert the structure to primitive form. This helps to
            reduce the number of fireworks in the workflow if the absorbing atoms is
            specified by its atomic symbol.
        feff_input_set (FeffDictSet)

    Returns:
        Workflow
    """
    if use_primitive:
        structure = structure.get_primitive_structure()

    # get the absorbing atom site index/indices
    ab_atom_indices = get_all_absorbing_atoms(absorbing_atom, structure)

    override_default_feff_params = {"user_tag_settings": user_tag_settings}

    # add firework for each absorbing atom site index
    fws = []
    for ab_idx in ab_atom_indices:
        fw_metadata = dict(metadata) if metadata else {}
        fw_metadata["absorbing_atom_index"] = ab_idx
        fw_name = "{}-{}-{}".format(spectrum_type, edge, ab_idx)
        fws.append(EELSFW(ab_idx, structure, spectrum_type, edge=edge, radius=radius,
                          beam_energy=beam_energy, beam_direction=beam_direction,
                          collection_angle=collection_angle, convergence_angle=convergence_angle,
                          user_eels_settings=user_eels_settings, feff_cmd=feff_cmd, db_file=db_file,
                          metadata=fw_metadata, name=fw_name, feff_input_set=feff_input_set,
                          override_default_feff_params=override_default_feff_params))

    wfname = "{}:{}:{} edge".format(structure.composition.reduced_formula,
                                    "{} spectroscopy".format(spectrum_type), edge)
    wf_metadata = dict(metadata) if metadata else {}
    wf_metadata["absorbing_atom_indices"] = list(ab_atom_indices)

    return Workflow(fws, name=wfname, metadata=wf_metadata)
