__version__ = '0.0.2'


from monty.json import MontyDecoder
from fireworks import Workflow, LaunchPad
from matmethods.vasp.vasp_powerups import decorate_write_name


def get_wf_from_spec_dict(structure, wfspec):
    fws = []
    common_params = wfspec.get("common_params", {})
    for d in wfspec["fireworks"]:
        modname, classname = d["fw"].rsplit(".", 1)
        mod = __import__(modname, globals(), locals(), [classname], 0)
        if hasattr(mod, classname):
            cls_ = getattr(mod, classname)
            params = {k: MontyDecoder().process_decoded(v) for k, v in d.get("params", {}).items()}
            params.update(common_params)
            if "parents" in params:
                params["parents"] = fws[params["parents"]]
            fws.append(cls_(structure, **params))
    return Workflow(fws, name=structure.composition.reduced_formula)


def add_to_lpad(workflow, decorate=False):
    """
    Add the workflow to the launchpad

    Args:
        workflow (Workflow): workflow for db insertion
        decorate (bool): If set an empty file with the name
            "FW--<fw.name>" will be written to the launch directory
    """
    lp = LaunchPad.auto_load()
    workflow = decorate_write_name(workflow) if decorate else workflow
    lp.add_wf(workflow)