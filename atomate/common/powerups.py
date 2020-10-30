"""
This module defines general powerups that can be used for all workflows
"""
from typing import Dict, Any, Optional

from atomate.utils.utils import get_fws_and_tasks
from fireworks.core.firework import Workflow

__author__ = "Janine George, Guido Petretto, Alex Ganose"
__email__ = "janine.george@uclouvain.be, guido.petretto@uclouvain.be"


def update_firetask(
    original_wf: Workflow,
    update_params: Dict[str, Any],
    task_name_constraint: str,
    fw_name_constraint: Optional[str] = None,
) -> Workflow:
    """
    General powerup for arbitrary updates to a Firetask.

    Args:
        original_wf: The original workflow.
        update_params: A dictionary of the keyword arguments to update.
        task_name_constraint: Only apply changes to the Firetasks where the
            Firetask class name contains this string.
        fw_name_constraint: Only apply changes to Fireworks where the Firework
            name contains this substring.

    Returns:
       Workflow
    """
    idx_list = get_fws_and_tasks(
        original_wf,
        fw_name_constraint=fw_name_constraint,
        task_name_constraint=task_name_constraint,
    )
    for idx_fw, idx_t in idx_list:
        original_wf.fws[idx_fw].tasks[idx_t].update(update_params)

    return original_wf


def set_queue_adapter(
    original_wf: Workflow,
    queueadapter: dict = None,

    fw_name_constraint: str = None,
    task_name_constraint: str = None,
) -> Workflow:
    """
    set _queueadapter spec of Fireworker(s) of a Workflow. It can be used to change the
    overall queueadapter during the run.

    Args:
        original_wf (Workflow): workflow that will be changed
        queueadapter (dict): dict to change _queueadapter
        fw_name_constraint (str): name of the Fireworks to be tagged (all if None is
            passed)
        task_name_constraint (str): name of the Firetasks to be tagged (e.g. None or
            'RunVasp')

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
