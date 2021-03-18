# coding: utf-8

"""
This module defines general powerups that can be used for all workflows
"""

from atomate.utils.utils import get_fws_and_tasks
from fireworks import Workflow, FileWriteTask
from fireworks.utilities.fw_utilities import get_slug


__author__ = "Janine George, Guido Petretto, Ryan Kingsbury"
__email__ = "janine.george@uclouvain.be, guido.petretto@uclouvain.be, RKingsbury@lbl.gov"


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
    Fireworks ending in "ToDb", e.g. VaspToDb, QChemToDb, etc., then the tags
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
        if original_wf.fws[idx_fw].tasks[idx_t].get("additional_fields"):
            if "tags" in original_wf.fws[idx_fw].tasks[idx_t]["additional_fields"]:
                original_wf.fws[idx_fw].tasks[idx_t]["additional_fields"][
                    "tags"
                ].extend(tags_list)
            else:
                original_wf.fws[idx_fw].tasks[idx_t]["additional_fields"][
                    "tags"
                ] = tags_list

    return original_wf


def add_namefile(original_wf, use_slug=True):
    """
    Every FireWork begins by writing an empty file with the name
    "FW--<fw.name>". This makes it easy to figure out what jobs are in what
    launcher directories, e.g. "ls -l launch*/FW--*" from within a "block" dir.

    Args:
        original_wf (Workflow)
        use_slug (bool): whether to replace whitespace-type chars with a slug

    Returns:
       Workflow
    """
    for idx, fw in enumerate(original_wf.fws):
        fname = "FW--{}".format(fw.name)
        if use_slug:
            fname = get_slug(fname)

        t = FileWriteTask(files_to_write=[{"filename": fname, "contents": ""}])
        original_wf.fws[idx].tasks.insert(0, t)
    return original_wf


def add_additional_fields_to_taskdocs(
    original_wf, update_dict=None, task_name_constraint="VaspToDb"
):
    """
    For all VaspToDbTasks in a given workflow, add information  to
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
