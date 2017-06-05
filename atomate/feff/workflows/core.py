# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines FEFF XAS(XANES/EXAFS) workflows.
"""

from fireworks import Workflow

from atomate.utils.utils import get_logger
from atomate.feff.fireworks.core import XASFW, EXAFSPathsFW, EELSFW
from atomate.feff.firetasks.write_inputs import get_feff_input_set_obj

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


def get_wf_xas(absorbing_atom, structure, feff_input_set="pymatgen.io.feff.sets.MPXANESSet",
               edge="K", radius=10.0, feff_cmd="feff", db_file=None, metadata=None,
               user_tag_settings=None, use_primitive=False):
    """
    Returns FEFF XANES/EXAFS spectroscopy workflow.

    Args:
        absorbing_atom (str/int): absorbing atom symbol or site index. If the symbol is given,
             then the returned workflow will have fireworks for each absorbing site with the
             same symbol.
        structure (Structure): input structure
        feff_input_set (str or FeffDictSet subclass): The inputset for setting params. If string
                then either the entire path to the class or spectrum type must be provided
                e.g. "pymatgen.io.feff.sets.MPXANESSet" or "XANES"
        edge (str): absorption edge. Example: K, L1, L2, L3
        radius (float): cluster radius in angstroms. Ignored for K space calculations
        feff_cmd (str): path to the feff binary
        db_file (str):  path to the db file.
        metadata (dict): meta data
        user_tag_settings (dict): override feff default tag settings
        use_primitive (bool): convert the structure to primitive form. This helps to
            reduce the number of fireworks in the workflow if the absorbing atom is
            specified by its atomic symbol.

    Returns:
        Workflow
    """
    if use_primitive:
        structure = structure.get_primitive_structure()

    # get the absorbing atom site index/indices
    ab_atom_indices = [absorbing_atom] if isinstance(absorbing_atom, int) else structure.indices_from_symbol(absorbing_atom)

    override_default_feff_params = {"user_tag_settings": user_tag_settings}

    spectrum_type = get_feff_input_set_obj(feff_input_set, ab_atom_indices[0], structure).__class__.__name__[2:-3]

    # add firework for each absorbing atom site index
    fws = []
    for ab_idx in ab_atom_indices:
        fw_metadata = dict(metadata) if metadata else {}
        fw_metadata["absorbing_atom_index"] = ab_idx
        fw_name = "{}-{}-{}".format(spectrum_type, edge, ab_idx)
        fws.append(XASFW(ab_idx, structure, edge=edge, radius=radius,
                         feff_input_set=feff_input_set, feff_cmd=feff_cmd, db_file=db_file,
                         metadata=fw_metadata, name=fw_name,
                         override_default_feff_params=override_default_feff_params))

    wf_metadata = dict(metadata) if metadata else {}
    wf_metadata["absorbing_atom_indices"] = list(ab_atom_indices)
    wfname = "{}:{}:{} edge".format(structure.composition.reduced_formula,
                                    "{} spectroscopy".format(spectrum_type), edge)

    return Workflow(fws, name=wfname, metadata=wf_metadata)


def get_wf_exafs_paths(absorbing_atom, structure, paths, degeneracies=None, edge="K", radius=10.0,
                       feff_input_set="pymatgen.io.feff.sets.MPEXAFSSet", feff_cmd="feff",
                       db_file=None, metadata=None, user_tag_settings=None, use_primitive=False,
                       labels=None, filepad_file=None):
    """
    Returns FEFF EXAFS spectroscopy workflow that generates the scattering amplitudes for the given
    list of scattering paths. The scattering amplitude output files(feffNNNN.dat files) are
    inserted to filepad(see fireworks.utilities.filepad.py) on completion.

    Args:
        absorbing_atom (str/int): absorbing atom symbol or site index. If the symbol is given,
             then the returned workflow will have fireworks for each absorbing site with the
             same symbol.
        structure (Structure): input structure
        paths (list): list of paths. path = list of site indices.
        degeneracies (list): list of path degeneracies.
        edge (str): absorption edge. Example: K, L1, L2, L3
        radius (float): cluster radius in angstroms. Ignored for K space calculations
        feff_input_set (str or FeffDictSet subclass): The inputset for setting params. If string
                then the entire path to the class must be provided
                e.g. "pymatgen.io.feff.sets.MPEXAFSSet"
        feff_cmd (str): path to the feff binary
        db_file (str):  path to the db file.
        metadata (dict): meta data
        user_tag_settings (dict): override feff default tag settings
        use_primitive (bool): convert the structure to primitive form. This helps to
            reduce the number of fireworks in the workflow if the absorbing atom is
            specified by its atomic symbol.
        labels ([str]): list of labels for the scattering amplitudes file contents inserted into
            filepad. Useful for fetching the data from filepad later.
        filepad_file (str): path to filepad connection settings file.

    Returns:
        Workflow
    """
    labels = labels or []
    wflow = get_wf_xas(absorbing_atom, structure, feff_input_set, edge, radius, feff_cmd,
                       db_file, metadata, user_tag_settings, use_primitive)
    paths_fw = EXAFSPathsFW(absorbing_atom, structure, paths, degeneracies=degeneracies, edge=edge,
                            radius=radius, name="EXAFS Paths", feff_input_set=feff_input_set,
                            feff_cmd=feff_cmd, labels=labels, filepad_file=filepad_file)
    # append the scattering paths firework to the regular EXAFS workflow.
    paths_wf = Workflow.from_Firework(paths_fw)
    wflow.append_wf(paths_wf, wflow.leaf_fw_ids)
    return wflow


def get_wf_eels(absorbing_atom, structure=None, feff_input_set="pymatgen.io.feff.sets.MPELNESSet",
                edge="K", radius=10., beam_energy=100, beam_direction=None, collection_angle=1,
                convergence_angle=1, user_eels_settings=None, user_tag_settings=None, feff_cmd="feff",
                db_file=None, metadata=None, use_primitive=False):
    """
    Returns FEFF ELNES/EXELFS spectroscopy workflow.

    Args:
        absorbing_atom (str): absorbing atom symbol
        structure (Structure): input structure. If None and mp_id is provided, the corresponding
            structure will be fetched from the Materials Project db.
        feff_input_set (str or FeffDictSet subclass): The inputset for setting params. If string
                then either the entire path to the class or spectrum type must be provided
                e.g. "pymatgen.io.feff.sets.MPELNESSet" or "ELNES"
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

    Returns:
        Workflow
    """
    if use_primitive:
        structure = structure.get_primitive_structure()

    # get the absorbing atom site index/indices
    ab_atom_indices = [absorbing_atom] if isinstance(absorbing_atom, int) else structure.indices_from_symbol(absorbing_atom)

    override_default_feff_params = {"user_tag_settings": user_tag_settings}

    spectrum_type = get_feff_input_set_obj(feff_input_set, ab_atom_indices[0], structure).__class__.__name__[2:-3]

    # add firework for each absorbing atom site index
    fws = []
    for ab_idx in ab_atom_indices:
        fw_metadata = dict(metadata) if metadata else {}
        fw_metadata["absorbing_atom_index"] = ab_idx
        fw_name = "{}-{}-{}".format(spectrum_type, edge, ab_idx)
        fws.append(EELSFW(ab_idx, structure, feff_input_set=feff_input_set, edge=edge, radius=radius,
                          beam_energy=beam_energy, beam_direction=beam_direction,
                          collection_angle=collection_angle, convergence_angle=convergence_angle,
                          user_eels_settings=user_eels_settings, feff_cmd=feff_cmd, db_file=db_file,
                          metadata=fw_metadata, name=fw_name,
                          override_default_feff_params=override_default_feff_params))

    wfname = "{}:{}:{} edge".format(structure.composition.reduced_formula,
                                    "{} spectroscopy".format(spectrum_type), edge)
    wf_metadata = dict(metadata) if metadata else {}
    wf_metadata["absorbing_atom_indices"] = list(ab_atom_indices)

    return Workflow(fws, name=wfname, metadata=wf_metadata)
