__version__ = '0.0.2'


from monty.json import MontyDecoder
from fireworks import Workflow, LaunchPad


def get_wf_from_spec_dict(structure, wfspec):
    fws = []
    for d in wfspec["fireworks"]:
        modname, classname = d["fw"].rsplit(".", 1)
        mod = __import__(modname, globals(), locals(), [classname], 0)
        if hasattr(mod, classname):
            cls_ = getattr(mod, classname)
            kwargs = {k: MontyDecoder().process_decoded(v) for k, v in d.get("params", {}).items()}
            if "parents" in kwargs:
                kwargs["parents"] = fws[kwargs["parents"]]
            fws.append(cls_(structure, **kwargs))
    return Workflow(fws, name=structure.composition.reduced_formula)