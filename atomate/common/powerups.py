from typing import Dict, Any, Optional

from fireworks import Workflow
from atomate.utils.utils import get_fws_and_tasks


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
