# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

import os

from atomate import get_wf_from_spec_dict
from monty.serialization import loadfn

__author__ = 'Anubhav Jain <ajain@lbl.gov>, Shyue Ping Ong <ongsp@eng.ucsd.edu>, Kiran Mathew <kmathew@lbl.gov>'


module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


def get_wf(structure, wf_filename, params=None, common_params=None, vis=None):
    """
    A generic function to load generic VASP library workflows, while
    overriding some of the parameters via the function arguments

    Args:
        structure: (Structure) structure to run
        wf_filename: filename in library subdir, e.g. "band_structure.yaml"
        params: (list of dicts) set params for each Firework; format is list
            that is same length as # of fws in the workflow
        common_params: (dict) set common params
        vis: (VaspInputSet) A VaspInputSet to use for the first FW

    Returns:
        A Workflow
    """
    d = loadfn(os.path.join(module_dir, "library", wf_filename))

    if params:
        if len(params) != len(d["fireworks"]):
            raise ValueError("The length of the params array must match the"
                             "length of the Fireworks array!")
        for idx, v in enumerate(params):
            if "params" not in d["fireworks"][idx]:
                d["fireworks"][idx]["params"] = {}
            d["fireworks"][idx]["params"].update(v)

    if common_params:
        if 'common_params' not in d:
            d["common_params"] = {}
        d["common_params"].update(common_params)

    if vis:
        if "params" not in d["fireworks"][0]:
            d["fireworks"][0]["params"] = {}
        d["fireworks"][0]["params"]["vasp_input_set"] = vis.as_dict()

    return get_wf_from_spec_dict(structure, d)