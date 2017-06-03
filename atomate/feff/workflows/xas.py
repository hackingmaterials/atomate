# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines FEFF XAS(XANES/EXAFS) workflows.
"""

from six import string_types

from fireworks import Workflow

from atomate.utils.utils import get_logger
from atomate.feff.fireworks.core import XASFW, EXAFSPathsFW

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


def get_wf_xas(absorbing_atom, structure, feff_input_set="XANES", edge="K", radius=10.0,
               feff_cmd="feff", db_file=None, metadata=None,
               user_tag_settings=None, use_primitive=False):
    """
    Returns FEFF XANES/EXAFS spectroscopy workflow.

    Args:
        absorbing_atom (str/int): absorbing atom symbol or site index. If the symbol is given,
             then the returned workflow will have fireworks for each absorbing site with the
             same symbol.
        feff_input_set (str or FeffDictSet subclass): The inputset for setting params. If string
                then either the entire path to the class or spectrum type must be provided
                e.g. "pymatgen.io.feff.sets.MPXANESSet" or "XANES"
        structure (Structure): input structure
        spectrum_type (str): XANES or EXAFS
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

    if isinstance(feff_input_set, string_types):
        if "pymatgen" in feff_input_set:
            _, classname = feff_input_set.strip().rsplit(".", 1)
        else:
            classname = "MP{}Set".format(feff_input_set)
            feff_input_set = "pymatgen.io.feff.sets.{}".format(classname)
    else:
        classname = feff_input_set.__class__.__name__
        
    spectrum_type = classname[2:-3]

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
                       feff_input_set=None, feff_cmd="feff", db_file=None, metadata=None,
                       user_tag_settings=None, use_primitive=False, labels=None, filepad_file=None):
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
        feff_input_set (FeffDictSet): the input set for the FEFF run
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
    feff_input_set = feff_input_set or "pymatgen.io.feff.sets.MPEXAFSSet"
    wflow = get_wf_xas(absorbing_atom, structure, feff_input_set, edge, radius, feff_cmd,
                       db_file, metadata, user_tag_settings, use_primitive)
    paths_fw = EXAFSPathsFW(absorbing_atom, structure, paths, degeneracies=degeneracies, edge=edge,
                            radius=radius, name="EXAFS Paths", feff_input_set=feff_input_set,
                            feff_cmd=feff_cmd, labels=labels, filepad_file=filepad_file)
    # append the scattering paths firework to the regular EXAFS workflow.
    paths_wf = Workflow.from_Firework(paths_fw)
    wflow.append_wf(paths_wf, wflow.leaf_fw_ids)
    return wflow
