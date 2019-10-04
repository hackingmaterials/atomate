# coding: utf-8


import os

from atomate.utils.utils import get_wf_from_spec_dict

from monty.serialization import loadfn

__author__ = 'Anubhav Jain, Shyue Ping Ong, Kiran Mathew'
__email__ = 'ajain@lbl.gov, ongsp@eng.ucsd.edu, kmathew@lbl.gov'


module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


def get_wf(structure, wf_filename, params=None, common_params=None, vis=None, wf_metadata=None):
    """
    Get a workflow given a structure and a name of file from the workflow library.

    Possible options for wf_filename are listed in: atomate.vasp.workflows.base.library and
    include band structure, dielectric constant, NEB, and more.

    You can also override some of the parameters via the function arguments.

    Args:
        structure: (Structure) structure to run
        wf_filename: filename in library subdir, e.g. "band_structure.yaml"
        params: (list of dicts) set params for each Firework; format is list
            that is same length as # of fws in the workflow
        common_params: (dict) set common params
        vis: (VaspInputSet) A VaspInputSet to use for the first FW
        wf_metadata: (dict) workflow metadata

    Returns:
        A Workflow
    """
    d = loadfn(os.path.join(module_dir, "library", wf_filename))

    if params:
        if len(params) != len(d["fireworks"]):
            raise ValueError("The length of the params array must match the length of the Fireworks array!")

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

    if wf_metadata:
        d["metadata"] = d.get("metadata", {})
        d["metadata"].update(wf_metadata)

    return get_wf_from_spec_dict(structure, d)
