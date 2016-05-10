# coding: utf-8
# Copyright (c) Materials Virtual Lab & HackingMaterials.
# Distributed under the terms of the BSD License.

from __future__ import division, unicode_literals, print_function

from monty.json import MontyDecoder
from fireworks import Workflow
import os

"""
Helper functions to load workflow from spec.
"""


def get_wf_from_spec_dict(structure, wfspec):
    """
    Load a WF from a structure and a spec dict. This allows simple
    custom workflows to be constructed quickly via a YAML file.

    Args:
        structure (Structure): An input structure object.
        wfspec (dict): A dict specifying workflow. A sample of the dict in
            YAML format for the usual MP workflow is given as follows:

            ```
            fireworks:
            - fw: matmethods.vasp.fireworks.core.OptimizeFW
            - fw: matmethods.vasp.fireworks.core.StaticFW
              params:
                parents: 0
            - fw: matmethods.vasp.fireworks.core.NonSCFUniformFW
              params:
                parents: 1
            - fw: matmethods.vasp.fireworks.core.NonSCFLineFW
              params:
                parents: 1
            common_params:
              db_file: db.json
              $vasp_cmd: $HOME/opt/vasp
            name: bandstructure
            ```

            The `fireworks` key is a list of fireworks. Each firework is
            specified via "fw": <explicit path>, and all parameters other than
            structure are specified via `params` which is a dict. `parents` is
            a special parameter, which provides the *indices* of the parents
            of that particular firework in the list. Any param starting with
            a $ will be expanded using environmental variables.

            `common_params` specify a common set of parameters that are
            passed to all fireworks, e.g., db settings.

            `name` is used to set the Workflow name (structure formula + name)

    Returns:
        Workflow
    """

    dec = MontyDecoder()

    def load_class(dotpath):
        modname, classname = dotpath.rsplit(".", 1)
        mod = __import__(modname, globals(), locals(), [classname], 0)
        return getattr(mod, classname)

    def process_params(d):
        decoded = {}
        for k, v in d.items():
            if k.startswith("$"):
                if isinstance(v, list):
                    v = [os.path.expandvars(i) for i in v]
                elif isinstance(v, dict):
                    v = {k2: os.path.expandvars(v2) for k2, v2 in v.items()}
                else:
                    v = os.path.expandvars(v)
            decoded[k.strip("$")] = dec.process_decoded(v)
        return decoded

    fws = []
    common_params = process_params(wfspec.get("common_params", {}))
    for d in wfspec["fireworks"]:
        cls_ = load_class(d["fw"])
        params = process_params(d.get("params", {}))
        params.update(common_params)
        if "parents" in params:
            if isinstance(params["parents"], int):
                params["parents"] = fws[params["parents"]]
            else:
                p = []
                for parent_idx in params["parents"]:
                    p.append(fws[parent_idx])
                params["parents"] = p
        fws.append(cls_(structure, **params))

    wfname = "{}:{}".format(structure.composition.reduced_formula, wfspec["name"]) if \
        wfspec.get("name") else structure.composition.reduced_formula
    return Workflow(fws, name=wfname)
