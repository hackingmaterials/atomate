"""
This module defines general powerups that can be used for all workflows
"""
from importlib import import_module
from typing import List

from fireworks import FileWriteTask, Workflow
from fireworks.utilities.fw_utilities import get_slug

from atomate.utils.utils import get_fws_and_tasks

__author__ = "Janine George, Guido Petretto, Ryan Kingsbury"
__email__ = (
    "janine.george@uclouvain.be, guido.petretto@uclouvain.be, RKingsbury@lbl.gov"
)


def add_priority(original_wf, root_priority, child_priority=None):
    """
    Adds priority to a workflow

    Args:
        original_wf (Workflow): original WF
        root_priority (int): priority of first (root) job(s)
        child_priority(int): priority of all child jobs. Defaults to
            root_priority

    Returns:
       Workflow: priority-decorated workflow
    """
    child_priority = child_priority or root_priority
    root_fw_ids = original_wf.root_fw_ids
    for fw in original_wf.fws:
        if fw.fw_id in root_fw_ids:
            fw.spec["_priority"] = root_priority
        else:
            fw.spec["_priority"] = child_priority
    return original_wf


def add_tags(original_wf, tags_list):
    """
    Adds tags to all Fireworks in the Workflow, WF metadata, as well as
    additional_fields for the Drone to track them later (e.g. tag all fireworks
    and tasks related to a specific research project).

    Tags are written to the "_spec" key of each Firework in the workflow and
    to the "metadata.tags" key of each Workflow. If the workflow contains any
    Firetasks ending in "ToDb", e.g. VaspToDb, QChemToDb, etc., then the tags
    are also passed as "additional_fields" to these tasks and included in the
    resulting task documents.

    Args:
        original_wf (Workflow)
        tags_list: list of tags parameters (list of strings)

    Returns:
       Workflow
    """

    # WF metadata
    if "tags" in original_wf.metadata:
        original_wf.metadata["tags"].extend(tags_list)
    else:
        original_wf.metadata["tags"] = tags_list

    # FW metadata
    for idx_fw in range(len(original_wf.fws)):
        if "tags" in original_wf.fws[idx_fw].spec:
            original_wf.fws[idx_fw].spec["tags"].extend(tags_list)
        else:
            original_wf.fws[idx_fw].spec["tags"] = tags_list

    # DB insertion tasks
    idxs = get_fws_and_tasks(original_wf, task_name_constraint="ToDb")
    for idx_fw, idx_t in idxs:
        if "additional_fields" in original_wf.fws[idx_fw].tasks[idx_t].optional_params:
            if "tags" in original_wf.fws[idx_fw].tasks[idx_t]["additional_fields"]:
                original_wf.fws[idx_fw].tasks[idx_t]["additional_fields"][
                    "tags"
                ].extend(tags_list)
            else:
                original_wf.fws[idx_fw].tasks[idx_t]["additional_fields"][
                    "tags"
                ] = tags_list

    return original_wf


def add_namefile(original_wf: Workflow, use_slug: bool = True) -> Workflow:
    """
    Every FireWork begins by writing an empty file with the name
    "FW--<fw.name>-<fw.fw_id>". This makes it easy to figure out what jobs are
    in what launcher directories, e.g. "ls -l launch*/FW--*" from within a
    "block" dir.

    Args:
        original_wf (Workflow)
        use_slug (bool): whether to replace whitespace-type chars with a slug.
            Defaults to True.

    Returns:
       Workflow
    """
    for idx, fw in enumerate(original_wf.fws):
        fname = f"FW--{fw.name}-{fw.fw_id}"
        if use_slug:
            fname = get_slug(fname)

        t = FileWriteTask(files_to_write=[{"filename": fname, "contents": ""}])
        original_wf.fws[idx].tasks.insert(0, t)
    return original_wf


def add_additional_fields_to_taskdocs(
    original_wf, update_dict=None, task_name_constraint="ToDb"
):
    """
    For all XXToDbTasks in a given workflow, add information  to
    "additional_fields" to be placed in the task doc.

    Args:
        original_wf (Workflow)
        update_dict (Dict): dictionary to add additional_fields
        task_name_constraint (str): name of the Firetasks to be modified.

    Returns:
       Workflow
    """
    idx_list = get_fws_and_tasks(original_wf, task_name_constraint=task_name_constraint)
    for idx_fw, idx_t in idx_list:
        original_wf.fws[idx_fw].tasks[idx_t]["additional_fields"].update(update_dict)
    return original_wf


def add_metadata(wf, meta_dict, fw_name_constraint=None):
    """
    Add a metadata dictionary to a Workflow and all its Fireworks. The dictionary
    is merged into the "metadata" key of the Workflow and into the "_spec" key of
    each Firework in the workflow.

    Can be used in combination with add_additional_fields_to_taskdocs to add the
    same set of key-value pairs to Workflows, Fireworks and Tasks collections.

    Args:
        wf (Workflow)
        meta_dict: dictionary of custom metadata

    Returns:
       Workflow
    """

    # add metadata to Workflow metadata
    wf.metadata.update(meta_dict)

    # add metadata to Firework metadata
    for fw in wf.fws:
        if fw_name_constraint is None or fw_name_constraint in fw.name:
            fw.spec.update(meta_dict)

    return wf


def preserve_fworker(original_wf, fw_name_constraint=None):
    """
    set _preserve_fworker spec of Fireworker(s) of a Workflow. Can be used to
    pin a workflow to the first fworker it is run with. Very useful when running
    on multiple machines that can't share files. fw_name_constraint can be used
    to only preserve fworker after a certain point where file passing becomes
    important

    Args:
        original_wf (Workflow):
        fw_name_constraint (str): name of the Fireworks to be tagged (all if
        None is passed)

    Returns:
        Workflow: modified workflow with specified Fireworkers tagged
    """
    idx_list = get_fws_and_tasks(original_wf, fw_name_constraint=fw_name_constraint)
    for idx_fw, idx_t in idx_list:
        original_wf.fws[idx_fw].spec["_preserve_fworker"] = True
    return original_wf


def set_execution_options(
    original_wf,
    fworker_name=None,
    category=None,
    fw_name_constraint=None,
    task_name_constraint=None,
):
    """
    set _fworker spec of Fireworker(s) of a Workflow. It can be used to specify
    a queue; e.g. run large-memory jobs on a separate queue.

    Args:
        original_wf (Workflow):
        fworker_name (str): user-defined tag to be added under fw.spec._fworker
            e.g. "large memory", "big", etc
        category (str): category of FWorker that should pul job
        fw_name_constraint (str): name of the Fireworks to be tagged (all if
            None is passed)
        task_name_constraint (str): name of the Firetasks to be tagged (e.g.
            None or 'RunVasp')

    Returns:
        Workflow: modified workflow with specified Fireworkers tagged
    """
    idx_list = get_fws_and_tasks(
        original_wf,
        fw_name_constraint=fw_name_constraint,
        task_name_constraint=task_name_constraint,
    )

    for idx_fw, idx_t in idx_list:
        if fworker_name:
            original_wf.fws[idx_fw].spec["_fworker"] = fworker_name
        if category:
            original_wf.fws[idx_fw].spec["_category"] = category
    return original_wf


def set_queue_adapter(
    original_wf: Workflow,
    queueadapter: dict = None,
    fw_name_constraint: str = None,
    task_name_constraint: str = None,
) -> Workflow:
    """
    set _queueadapter spec of Fireworker(s) of a Workflow. It can be used to change the overall queueadapter during the run.

    Args:
        original_wf (Workflow): workflow that will be changed
        queueadapter (dict): dict to change _queueadapter
        fw_name_constraint (str): name of the Fireworks to be tagged (all if None is passed)
        task_name_constraint (str): name of the Firetasks to be tagged (e.g. None or 'RunVasp')

    Returns:
        Workflow: modified workflow with specified Fireworkers tagged
    """
    for idx_fw, idx_t in get_fws_and_tasks(
        original_wf,
        fw_name_constraint=fw_name_constraint,
        task_name_constraint=task_name_constraint,
    ):
        q = original_wf.fws[idx_fw].spec.get("_queueadapter", {})
        q.update(queueadapter)
        original_wf.fws[idx_fw].spec["_queueadapter"] = q

    return original_wf


def powerup_by_kwargs(
    original_wf: Workflow,
    powerup_dicts: List[dict],
):
    """
    apply powerups in the form using a list of dictionaries
    [
        {"powerup_name" : powerup_function1, "kwargs": {parameter1 : value1, parameter2: value2}},
        {"powerup_name" : powerup_function2, "kwargs": {parameter1 : value1, parameter2: value2}},
    ]

    As an example:
        power_up_by_kwargs([
            {"powerup_name" : "add_additional_fields_to_taskdocs",
             "kwargs: {"update_dict" : {"foo" : "bar"}}}
             ]
        )

    Args:
        original_wf: workflow that will be changed
        powerup_dicts: dictionary containing the powerup_name and kwarg.
            if "." is present in the name it will be imported as a full path
            if not we will use standard atomate modules where the powerups are kept

    """
    # a list of possible powerups in atomate (most specific first)
    powerup_modules = [
        "atomate.vasp.powerups",
        "atomate.qchem.powerups",
        "atomate.common.powerups",
    ]

    for pd in powerup_dicts:
        name = pd["powerup_name"]
        kwargs = pd["kwargs"]
        found = False
        if "." in name:
            module_name, method_name = name.rsplit(".", 1)
            module = import_module(module_name)
            powerup = getattr(module, method_name)
            original_wf = powerup(original_wf, **kwargs)
            found = True
        else:
            for module_name in powerup_modules:
                try:
                    module = import_module(module_name)
                    powerup = getattr(module, name)
                    original_wf = powerup(original_wf, **kwargs)
                    found = True
                    break
                except Exception:
                    pass
        if not found:
            raise RuntimeError(f"Could not find powerup {name}.")
    return original_wf
