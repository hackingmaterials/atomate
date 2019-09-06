# coding: utf-8


import logging
import os
import sys
import socket
from random import randint
from time import time

from pymongo import MongoClient
from monty.json import MontyDecoder
from monty.serialization import loadfn
from pymatgen import Composition

from fireworks import Workflow
from pymatgen.alchemy.materials import TransformedStructure

__author__ = 'Anubhav Jain, Kiran Mathew'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov'


def env_chk(val, fw_spec, strict=True, default=None):
    """
    env_chk() is a way to set different values for a property depending
    on the worker machine. For example, you might have slightly different
    executable names or scratch directories on different machines.

    env_chk() works using the principles of the FWorker env in FireWorks.

    This helper method translates string "val" that looks like this:
    ">>ENV_KEY<<"
    to the contents of:
    fw_spec["_fw_env"][ENV_KEY]

    Otherwise, the string "val" is interpreted literally and passed-through as is.

    The fw_spec["_fw_env"] is in turn set by the FWorker. For more details,
    see: https://materialsproject.github.io/fireworks/worker_tutorial.html

    Since the fw_env can be set differently for each FireWorker, one can
    use this method to translate a single "val" into multiple possibilities,
    thus achieving different behavior on different machines.

    Args:
        val: any value, with ">><<" notation reserved for special env lookup values
        fw_spec: (dict) fw_spec where one can find the _fw_env keys
        strict (bool): if True, errors if env format (>><<) specified but cannot be found in fw_spec
        default: if val is None or env cannot be found in non-strict mode,
                 return default
    """
    if val is None:
        return default

    if isinstance(val, str) and val.startswith(">>") and val.endswith("<<"):
        if strict:
            return fw_spec['_fw_env'][val[2:-2]]
        return fw_spec.get('_fw_env', {}).get(val[2:-2], default)
    return val


def get_mongolike(d, key):
    """
    Retrieve a dict value using dot-notation like "a.b.c" from dict {"a":{"b":{"c": 3}}}
    Args:
        d (dict): the dictionary to search
        key (str): the key we want to retrieve with dot notation, e.g., "a.b.c"

    Returns:
        value from desired dict (whatever is stored at the desired key)

    """
    lead_key = key.split(".", 1)[0]
    try:
        lead_key = int(lead_key)  # for searching array data
    except:
        pass

    if "." in key:
        remainder = key.split(".", 1)[1]
        return get_mongolike(d[lead_key], remainder)
    return d[lead_key]


def recursive_get_result(d, result):
    """
    Function that gets designated keys or values of d
    (i. e. those that start with "d>>" or "a>>") from
    the corresponding entry in result_dict, similar to
    FireWorks recursive_deserialize.

    Note that the plain ">>" notation will get a key from
    the result.as_dict() object and may use MongoDB
    dot notation, while "a>>" will get an attribute
    of the object.

    Examples:

    Getting a dict key from a VaspRun instance:
        recursive_get_result({"stress":">>output.ionic_steps.-1.stress"}, vasprun)
        --> {"stress":[[0.2, 0, 0], [0, 0.3, 0], [0, 0, 0.3]]}

    Getting an **attribute** from a vasprun:
        recursive_get_result({"epsilon":"a>>epsilon_static", vasprun}
        --> {"epsilon":-3.4}
    """
    if isinstance(d, str) and d[:2] == ">>":
        if hasattr(result, "as_dict"):
            result = result.as_dict()
        return get_mongolike(result, d[2:])

    elif isinstance(d, str) and d[:3] == "a>>":
        attribute = getattr(result, d[3:])
        if callable(attribute):
            attribute = attribute()
        return attribute

    elif isinstance(d, dict):
        return {k: recursive_get_result(v, result) for k, v in d.items()}

    elif isinstance(d, (list, tuple)):
        return [recursive_get_result(i, result) for i in d]

    else:
        return d


def get_logger(name, level=logging.DEBUG, log_format='%(asctime)s %(levelname)s %(name)s %(message)s', stream=sys.stdout):
    logger = logging.getLogger(name)
    logger.setLevel(level)
    formatter = logging.Formatter(log_format)
    sh = logging.StreamHandler(stream=stream)
    sh.setFormatter(formatter)
    logger.addHandler(sh)
    return logger


def get_meta_from_structure(structure):
    if isinstance(structure, TransformedStructure):
        structure = structure.final_structure

    comp = structure.composition
    elsyms = sorted(set([e.symbol for e in comp.elements]))
    meta = {'nsites': structure.num_sites,
            'elements': elsyms,
            'nelements': len(elsyms),
            'formula': comp.formula,
            'formula_pretty': comp.reduced_formula,
            'formula_reduced_abc': Composition(comp.reduced_formula)
            .alphabetical_formula,
            'formula_anonymous': comp.anonymized_formula,
            'chemsys': '-'.join(elsyms),
            'is_ordered': structure.is_ordered,
            'is_valid': structure.is_valid()}
    return meta


def get_fws_and_tasks(workflow, fw_name_constraint=None, task_name_constraint=None):
    """
    Helper method: given a workflow, returns back the fw_ids and task_ids that match name
    constraints. Used in developing multiple powerups.

    Args:
        workflow (Workflow): Workflow
        fw_name_constraint (str): a constraint on the FW name
        task_name_constraint (str): a constraint on the task name

    Returns:
       a list of tuples of the form (fw_id, task_id) of the RunVasp-type tasks
    """
    fws_and_tasks = []
    for idx_fw, fw in enumerate(workflow.fws):
        if fw_name_constraint is None or fw_name_constraint in fw.name:
            for idx_t, t in enumerate(fw.tasks):
                if task_name_constraint is None or task_name_constraint in str(t):
                    fws_and_tasks.append((idx_fw, idx_t))
    return fws_and_tasks


# TODO: @computron - move this somewhere else, maybe dedicated serialization package - @computron
# TODO: @computron - also review this code for clarity - @computron
def get_wf_from_spec_dict(structure, wfspec, common_param_updates=None):
    """
    Load a WF from a structure and a spec dict. This allows simple
    custom workflows to be constructed quickly via a YAML file.

    Args:
        structure (Structure): An input structure object.
        wfspec (dict): A dict specifying workflow. A sample of the dict in
            YAML format for the usual MP workflow is given as follows:

            ```
            fireworks:
            - fw: atomate.vasp.fireworks.core.OptimizeFW
            - fw: atomate.vasp.fireworks.core.StaticFW
              params:
                parents: 0
            - fw: atomate.vasp.fireworks.core.NonSCFUniformFW
              params:
                parents: 1
            - fw: atomate.vasp.fireworks.core.NonSCFLineFW
              params:
                parents: 1
            common_params:
              db_file: db.json
              $vasp_cmd: $HOME/opt/vasp
            name: bandstructure
            metadata:
                tag: testing_workflow
            ```

            The `fireworks` key is a list of Fireworks; it is expected that
            all such Fireworks have "structure" as the first argument and
            other optional arguments following that. Each Firework is specified
            via "fw": <explicit path>.

            You can pass arguments into the constructor using the special
            keyword `params`, which is a dict. Any param starting with a $ will
            be expanded using environment variables.If multiple fireworks share
            the same `params`, you can use `common_params` to specify a common
            set of arguments that are passed to all fireworks. Local params
            take precedent over global params.

            Another special keyword is `parents`, which provides
            the *indices* of the parents of that particular Firework in the
            list. This allows you to link the Fireworks into a logical
            workflow.

            Finally, `name` is used to set the Workflow name
            (structure formula + name) which can be helpful in record keeping.
        common_param_updates (dict): A dict specifying any user-specified updates to common_params

    Returns:
        Workflow
    """

    dec = MontyDecoder()

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
    if common_param_updates:
        common_params.update(common_param_updates)
    for d in wfspec["fireworks"]:
        modname, classname = d["fw"].rsplit(".", 1)
        cls_ = load_class(modname, classname)
        params = process_params(d.get("params", {}))
        for k in common_params:
            if k not in params:  # common params don't override local params
                params[k] = common_params[k]
        if "parents" in params:
            if isinstance(params["parents"], int):
                params["parents"] = fws[params["parents"]]
            else:
                p = []
                for parent_idx in params["parents"]:
                    p.append(fws[parent_idx])
                params["parents"] = p
        fws.append(cls_(structure=structure, **params))

    wfname = "{}:{}".format(structure.composition.reduced_formula, wfspec["name"]) if \
        wfspec.get("name") else structure.composition.reduced_formula

    return Workflow(fws, name=wfname, metadata=wfspec.get("metadata"))


def load_class(modulepath, classname):
    """
    Load and return the class from the given module.

    Args:
        modulepath (str): dotted path to the module. eg: "pymatgen.io.vasp.sets"
        classname (str): name of the class to be loaded.

    Returns:
        class
    """
    mod = __import__(modulepath, globals(), locals(), [classname], 0)
    return getattr(mod, classname)


def recursive_update(d, u):
    """
    Recursive updates d with values from u
    Args:
        d (dict): dict to update
        u (dict): updates to propogate
    """

    for k, v in u.items():
        if k in d:
            if isinstance(v, dict) and isinstance(d[k], dict):
                recursive_update(d[k], v)
            else:
                d[k] = v
        else:
            d[k] = v


def get_a_unique_id():
    ts = "{:.4f}".format(time())
    ts += str(randint(0, 9999)).zfill(4)
    return ts


def get_uri(dir_name):
    """
    Returns the URI path for a directory. This allows files hosted on
    different file servers to have distinct locations.
    Args:
        dir_name:
            A directory name.
    Returns:
        Full URI path, e.g., fileserver.host.com:/full/path/of/dir_name.
    """
    fullpath = os.path.abspath(dir_name)
    try:
        hostname = socket.gethostbyaddr(socket.gethostname())[0]
    except:
        hostname = socket.gethostname()
    return "{}:{}".format(hostname, fullpath)


def get_database(config_file=None, settings=None, admin=False, **kwargs):
    d = loadfn(config_file) if settings is None else settings

    try:
        user = d["admin_user"] if admin else d["readonly_user"]
        passwd = d["admin_password"] if admin else d["readonly_password"]
    except (KeyError, TypeError, ValueError):
        logger.warning("No {admin,readonly}_user/password found in config. file, "
            "accessing DB without authentication")
        user = None
        passwd = None

    if "authsource" in d and "authsource" not in kwargs:
        kwargs["authsource"] = d["authsource"]

    conn = MongoClient(host=d["host"], port=d["port"], username=user,
                       password=passwd, **kwargs)
    db = conn[d["database"]]

    return db


logger = get_logger(__name__)
