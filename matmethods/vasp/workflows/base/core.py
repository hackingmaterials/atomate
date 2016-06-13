import os

from matmethods import get_wf_from_spec_dict
from monty.serialization import loadfn
from pymatgen.io.vasp.sets import MPRelaxSet

__author__ = 'Anubhav Jain <ajain@lbl.gov>, ' \
             'Shyue Ping Ong <ongsp@eng.ucsd.edu>, ' \
             'Kiran Mathew <kmathew@lbl.gov>'


module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


def get_wf(structure, wf_filename, vasp_input_set=None, vasp_cmd="vasp",
           db_file=None, params=None, common_params=None):

    """
    A generic function to load most library workflows (see library subdir)
    Args:
        structure: (Structure) structure to compute
        wf_filename: filename in library subdir, e.g. "band_structure.yaml"
        vasp_input_set: VaspInputSet for first job (usually chains)
        vasp_cmd: (str) vasp command (use >><< notation for env_chk)
        db_file: (str) db_file (use >><< notation for env_chk)
        params: (list of dicts) set params for each Firework; format is list
            that is same length as # of fws in the workflow
        common_params: (dict) set common params

    Returns:
        A Workflow

    """
    d = loadfn(os.path.join(module_dir, "library", wf_filename))
    v = vasp_input_set or MPRelaxSet(structure, force_gamma=True)
    if "params" not in d["fireworks"][0]:
        d["fireworks"][0]["params"] = {}
    d["fireworks"][0]["params"]["vasp_input_set"] = v.as_dict()
    d["common_params"] = {"vasp_cmd": vasp_cmd, "db_file": db_file}

    if params:
        if len(params) != len(d["fireworks"]):
            raise ValueError("The length of the params array must match the"
                             "length of the Fireworks array!")
        for idx, v in enumerate(params):
            d["fireworks"][idx].update(v)

    if common_params:
        d["common_params"].update(common_params)

    return get_wf_from_spec_dict(structure, d)