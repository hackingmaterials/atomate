# coding: utf-8
# Copyright (c) Materials Virtual Lab & HackingMaterials.
# Distributed under the terms of the BSD License.

from __future__ import division, unicode_literals, print_function

"""
Various helper functions to load and add workflows from YAML specs.
"""

from monty.json import MontyDecoder
from fireworks import Workflow


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
              vasp_cmd: /opt/vasp
            ```

            The `fireworks` key is a list of fireworks. Each firework is
            specified via "fw": <explicit path>, and all parameters other than
            structure are specified via `params` which is a dict. `parents` is
            a special parameter, which provides the *indices* of the parents
            of that particular firework in the list.

            `common_params` specify a common set of parameters that are
            passed to all fireworks, e.g., db settings.

    Returns:
        Workflow
    """
    fws = []
    common_params = wfspec.get("common_params", {})
    for d in wfspec["fireworks"]:
        modname, classname = d["fw"].rsplit(".", 1)
        mod = __import__(modname, globals(), locals(), [classname], 0)
        cls_ = getattr(mod, classname)
        params = {k: MontyDecoder().process_decoded(v) for k, v in d.get("params", {}).items()}
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

    return Workflow(fws, name=structure.composition.reduced_formula)
