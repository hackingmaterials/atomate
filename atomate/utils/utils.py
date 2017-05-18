# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import logging
import os
import sys
from copy import deepcopy

import six

from fireworks import Workflow, Firework
from monty.json import MontyDecoder
from pymatgen import Composition

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
    see: https://pythonhosted.org/FireWorks/worker_tutorial.html

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

    if isinstance(val, six.string_types) and val.startswith(">>") and val.endswith("<<"):
        if strict:
            return fw_spec['_fw_env'][val[2:-2]]
        return fw_spec.get('_fw_env', {}).get(val[2:-2], default)
    return val


def get_mongolike(d, key):
    """
    Grab a dict value using dot-notation like "a.b.c" from dict {"a":{"b":{"c": 3}}}
    Args:
        d (dict): the dictionary to search
        key (str): the key we want to grab with dot notation, e.g., "a.b.c" 

    Returns:
        value from desired dict (whatever is stored at the desired key)

    """
    if "." in key:
        i, j = key.split(".", 1)
        try:
            i = int(i)  # for searching array data
        except:
            pass
        return get_mongolike(d[i], j)
    return d[key]


def get_logger(name, level=logging.DEBUG, format='%(asctime)s %(levelname)s %(name)s %(message)s',
               stream=sys.stdout):
    logger = logging.getLogger(name)
    logger.setLevel(level)
    formatter = logging.Formatter(format)
    sh = logging.StreamHandler(stream=stream)
    sh.setFormatter(formatter)
    logger.addHandler(sh)
    return logger


def get_meta_from_structure(structure):
    comp = structure.composition
    elsyms = sorted(set([e.symbol for e in comp.elements]))
    meta = {'nsites': len(structure),
            'elements': elsyms,
            'nelements': len(elsyms),
            'formula': comp.formula,
            'formula_reduced': comp.reduced_formula,
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
        fws.append(cls_(structure, **params))

    wfname = "{}:{}".format(structure.composition.reduced_formula, wfspec["name"]) if \
        wfspec.get("name") else structure.composition.reduced_formula

    return Workflow(fws, name=wfname, metadata=wfspec.get("metadata"))


# TODO: @matk86 - please remove this pointless method. Write tighter code rather than this silly
# "auto-converting without thinking about it" mess. Even if you wanted to do this auto-conversion,
# just modify Workflow.append_wf() in FireWorks rather than add this thing. -computron
def append_fw_wf(orig_wf, fw_wf):
    """
    Add the given firework or workflow to the end of the provided workflow. If there are multiple
    leaf nodes the newly added firework/workflow will depend on all of them.

    Args:
        orig_wf (Workflow): The original workflow object.
        fw_wf (Firework/Workflow): The firework or workflow object to be appended to orig_wf.
    """
    new_wf = fw_wf
    if isinstance(fw_wf, Firework):
        new_wf = Workflow.from_Firework(new_wf)
    orig_wf.append_wf(new_wf, orig_wf.leaf_fw_ids)


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


# TODO: @matk86 - put this in FireWorks instead? w/ already existing unit tests. -computron
def remove_fws(orig_wf, fw_ids):
    """
    Remove the fireworks corresponding to the input firework ids and update the workflow i.e the
    parents of the removed fireworks become the parents of the children fireworks (only if the
    children dont have any other parents).

    Args:
        orig_wf (Workflow): The original workflow object.
        fw_ids (list): list of fw ids to remove.

    Returns:
        Workflow : the new updated workflow.
    """
    # not working with the copies causes spurious behavior
    wf_dict = deepcopy(orig_wf.as_dict())
    orig_parent_links = deepcopy(orig_wf.links.parent_links)
    fws = wf_dict["fws"]

    # update the links dict: remove fw_ids and link their parents to their children (if they don't
    # have any other parents).
    for fid in fw_ids:
        children = wf_dict["links"].pop(str(fid))
        # root node --> no parents
        try:
            parents = orig_parent_links[int(fid)]
        except KeyError:
            parents = []
        # remove the firework from their parent links and re-link their parents to the children.
        for p in parents:
            wf_dict["links"][str(p)].remove(fid)
            # adopt the children
            for c in children:
                # adopt only if the child doesn't have any other parents.
                if len(orig_parent_links[int(c)]) == 1:
                    wf_dict["links"][str(p)].append(c)

    # update the list of fireworks.
    wf_dict["fws"] = [f for f in fws if f["fw_id"] not in fw_ids]

    new_wf = Workflow.from_dict(wf_dict)

    return new_wf
