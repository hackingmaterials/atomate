# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

from fireworks import Workflow, FileWriteTask
from fireworks.core.firework import Tracker
from fireworks.utilities.fw_utilities import get_slug
from matmethods.vasp.firetasks.glue_tasks import CheckStability

from matmethods.vasp.firetasks.run_calc import RunVaspCustodian, RunVaspDirect
from matmethods.vasp.firetasks.write_inputs import ModifyIncar
from matmethods.vasp.tests.vasp_fake import fake_dirs, RunVaspFake

__author__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'


def get_fws_and_tasks(workflow, fw_name_constraint=None, task_name_constraint=None):
    """
    Helper method: given a workflow, returns back the fw_ids and task_ids that match constraints

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


def add_priority(original_wf, root_priority, child_priority=None):
    """
    Adds priority to a workflow

    Args:
        original_wf (Workflow): original WF
        root_priority (int): priority of first (root) job(s)
        child_priority(int): priority of all child jobs. Defaults to
                            root_priority

    Returns:
       (Workflow) priority-decorated workflow
    """

    child_priority = child_priority or root_priority
    root_fw_ids = original_wf.root_fw_ids
    for fw in original_wf.fws:
        if fw.fw_id in root_fw_ids:
            fw.spec["_priority"] = root_priority
        else:
            fw.spec["_priority"] = child_priority
    return original_wf


def remove_custodian(original_wf, fw_name_constraint=None):
    """
    Replaces all tasks with "RunVasp*" (e.g. RunVaspCustodian) to be
    RunVaspDirect.

    Args:
        original_wf (Workflow): original workflow
        fw_name_constraint (str): Only apply changes to FWs where fw_name
            contains this substring.
    """

    wf_dict = original_wf.to_dict()
    vasp_fws_and_tasks = get_fws_and_tasks(original_wf, fw_name_constraint=fw_name_constraint,
                                           task_name_constraint="RunVasp")

    for idx_fw, idx_t in vasp_fws_and_tasks:
        vasp_cmd = wf_dict["fws"][idx_fw]["spec"]["_tasks"][idx_t]["vasp_cmd"]
        wf_dict["fws"][idx_fw]["spec"]["_tasks"][idx_t] = \
            RunVaspDirect(vasp_cmd=vasp_cmd).to_dict()

    return Workflow.from_dict(wf_dict)


def use_custodian(original_wf, fw_name_constraint=None, custodian_params=None):
    """
    Replaces all tasks with "RunVasp*" (e.g. RunVaspDirect) to be
    RunVaspCustodian. Thus, this powerup adds error correction into VASP
    runs if not originally present and/or modifies the correction behavior.

    Args:
        original_wf (Workflow): original workflow
        fw_name_constraint (str): Only apply changes to FWs where fw_name contains this substring.
            For example, use custodian only for certain runs, or set job_type to
            "double_relaxation_run" only for structure optimization run, or set different
            handler_group for different runs.
        custodian_params (dict): A dict of parameters for RunVaspCustodian. e.g., use it to set
            a "scratch_dir" or "handler_group".
    """

    custodian_params = custodian_params if custodian_params else {}
    wf_dict = original_wf.to_dict()
    vasp_fws_and_tasks = get_fws_and_tasks(original_wf, fw_name_constraint=fw_name_constraint,
                                           task_name_constraint="RunVasp")

    for idx_fw, idx_t in vasp_fws_and_tasks:
        if "vasp_cmd" not in custodian_params:
            custodian_params["vasp_cmd"] = \
                wf_dict["fws"][idx_fw]["spec"]["_tasks"][idx_t]["vasp_cmd"]

        wf_dict["fws"][idx_fw]["spec"]["_tasks"][idx_t] = \
            RunVaspCustodian(**custodian_params).to_dict()

    return Workflow.from_dict(wf_dict)


def use_fake_vasp(original_wf):
    """
    Replaces all tasks with "RunVasp" (e.g. RunVaspDirect) to be
    RunVaspFake. Thus, we do not actually run VASP but copy
    pre-determined inputs and outputs.

    Args:
        original_wf (Workflow)
    """
    wf_dict = original_wf.to_dict()
    for idx_fw, fw in enumerate(original_wf.fws):
        for job_type in fake_dirs.keys():
            if job_type in fw.name:
                for idx_t, t in enumerate(fw.tasks):
                    if "RunVasp" in str(t):
                        wf_dict["fws"][idx_fw]["spec"]["_tasks"][idx_t] = \
                            RunVaspFake(fake_dir=fake_dirs[job_type]).to_dict()

    return Workflow.from_dict(wf_dict)


def add_namefile(original_wf, use_slug=True):
    """
    Every FireWork begins by writing an empty file with the name
    "FW--<fw.name>". This makes it easy to figure out what jobs are in what
    launcher directories, e.g. "ls -l launch*/FW--*" from within a "block" dir.

    Args:
        original_wf (Workflow)
        use_slug (bool): whether to replace whitespace-type chars with a slug
    """
    wf_dict = original_wf.to_dict()

    for idx, fw in enumerate(wf_dict["fws"]):
        fname = "FW--{}".format(fw["name"])
        if use_slug:
            fname = get_slug(fname)
        wf_dict["fws"][idx]["spec"]["_tasks"].insert(0, FileWriteTask(
            files_to_write=[{"filename": fname, "contents": ""}]).to_dict())

    return Workflow.from_dict(wf_dict)


def add_trackers(original_wf):
    """
    Every FireWork that runs VASP also tracks the OUTCAR and OSZICAR using FWS Trackers.

    Args:
        original_wf (Workflow)

    """
    tracker1 = Tracker('OUTCAR', nlines=25, allow_zipped=True)
    tracker2 = Tracker('OSZICAR', nlines=25, allow_zipped=True)

    wf_dict = original_wf.to_dict()

    for idx_fw, idx_t in get_fws_and_tasks(original_wf, task_name_constraint="RunVasp"):
        if "_trackers" in wf_dict["fws"][idx_fw]["spec"]:
            wf_dict["fws"][idx_fw]["spec"]["_trackers"].extend([tracker1, tracker2])
        else:
            wf_dict["fws"][idx_fw]["spec"]["_trackers"] = [tracker1, tracker2]

    return Workflow.from_dict(wf_dict)


def add_modify_incar(original_wf, modify_incar_params=None, fw_name_constraint=None):
    """
    Every FireWork that runs VASP has a ModifyIncar task just beforehand. For example, allows
    you to modify the INCAR based on the Worker using env_chk or using hard-coded changes.

    Args:
        original_wf (Workflow)
        modify_incar_params (dict) - dict of parameters for ModifyIncar.
        fw_name_constraint (str) - Only apply changes to FWs where fw_name contains this substring.

    """

    modify_incar_params = modify_incar_params or {"incar_update": ">>incar_update<<"}
    for idx_fw, idx_t in get_fws_and_tasks(original_wf, fw_name_constraint=fw_name_constraint,
                                           task_name_constraint="RunVasp"):
        original_wf.fws[idx_fw].spec["_tasks"].insert(idx_t, ModifyIncar(**modify_incar_params).
                                                      to_dict())

    return original_wf


def add_stability_check(original_wf, check_stability_params=None, fw_name_constraint=None):
    """
    Every FireWork that runs VASP has a CheckStability task afterward. This
    allows defusing jobs that are not stable. In practice, you might want
    to set the fw_name_constraint so that the stability is only checked at the
    beginning of the workflow

    Args:
        original_wf (Workflow)
        check_stability_params (dict): a **kwargs** style dict of params
        fw_name_constraint (str) - Only apply changes to FWs where fw_name contains this substring.
    """

    check_stability_params = check_stability_params or {}
    for idx_fw, idx_t in get_fws_and_tasks(original_wf, fw_name_constraint=fw_name_constraint,
                                           task_name_constraint="RunVasp"):
        original_wf.fws[idx_fw].spec["_tasks"].append(CheckStability(**check_stability_params).to_dict())

    return original_wf


def add_modify_incar_envchk(original_wf, fw_name_constraint=None):
    """
    If you set the "incar_update" parameter in the Worker env, the INCAR will update this
    parameter for all matching VASP runs

    Args:
        original_wf (Workflow)
        fw_name_constraint (str) - Only apply changes to FWs where fw_name contains this substring.
    """
    return add_modify_incar(original_wf, {"incar_update": ">>incar_update<<"},
                            fw_name_constraint=fw_name_constraint)


def add_small_gap_multiply(original_wf, gap_cutoff, density_multiplier, fw_name_constraint=None):
    """
    In all FWs with specified name constraints, add a 'small_gap_multiply' parameter that
    multiplies the k-mesh density of compounds with gap < gap_cutoff by density multiplier.
    Note that this powerup only works on FireWorks with the appropriate WriteVasp* tasks that
    accept the small_gap_multiply argument...

    :param original_wf:
    :param gap_cutoff:
    :param density_multiplier:
    :param fw_name_constraint:
    """

    wf_dict = original_wf.to_dict()
    for idx_fw, idx_t in get_fws_and_tasks(original_wf, fw_name_constraint=fw_name_constraint,
                                           task_name_constraint="WriteVasp"):
        wf_dict["fws"][idx_fw]["spec"]["_tasks"][idx_t]["small_gap_multiply"] = \
            [gap_cutoff, density_multiplier]

    return Workflow.from_dict(wf_dict)


def use_scratch_dir(original_wf, scratch_dir):
    """
    For all RunVaspCustodian tasks, add the desired scratch dir.

    :param original_wf:
    :param scratch_dir: The scratch dir to use. Supports env_chk
    """

    wf_dict = original_wf.to_dict()
    for idx_fw, idx_t in get_fws_and_tasks(original_wf, task_name_constraint="RunVaspCustodian"):
        wf_dict["fws"][idx_fw]["spec"]["_tasks"][idx_t]["scratch_dir"] = scratch_dir

    return Workflow.from_dict(wf_dict)

def modify_task_docs(original_wf, update_dict = None):
    """
    For all VaspToDbTasks in a given workflow, add information 
    to "additional_fields" to be placed in the task doc.
    
    Args:
        original_wf (Workflow)
        update_dict (Dict): dictionary to add additional_fields
    """
    wf_dict = original_wf.to_dict()
    for idx_fw, idx_t in get_fws_and_tasks(original_wf, 
                                           task_name_constraint="VaspToDbTask"):
        wf_dict["fws"][idx_fw]["spec"]["_tasks"][idx_t]["additional_fields"].update(
            update_dict)
    return Workflow.from_dict(wf_dict)

