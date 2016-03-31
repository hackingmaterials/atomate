# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

from fireworks import Workflow, FileWriteTask
from fireworks.core.firework import Tracker
from fireworks.utilities.fw_utilities import get_slug

from matmethods.vasp.firetasks.run_calc import RunVaspCustodian
from matmethods.vasp.tests.vasp_fake import fake_dirs, RunVaspFake

__author__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'


def is_runvasp_task(task):
    """
    Given a FireTask, tells you if it runs VASP.

    Args:
        task (FireTask): FireTask

    Returns:
       True/False
    """
    if "RunVasp" in str(task):
        return True

    return False


def get_runvasp_fws_and_tasks(workflow):
    """
    Given a workflow, returns back the fw_ids and task_ids of RunVasp-type tasks

    Args:
        workflow (Workflow): Workflow

    Returns:
       a list of tuples of the form (fw_id, task_id) of the RunVasp-type tasks
    """
    fws_and_tasks = []

    for idx_fw, fw in enumerate(workflow.fws):
        for job_type in fake_dirs.keys():
            if job_type in fw.name:
                for idx_t, t in enumerate(fw.tasks):
                    if is_runvasp_task(t):
                        fws_and_tasks.append((idx_fw, idx_t))

    return fws_and_tasks


def decorate_priority(original_wf, root_priority, child_priority=None):
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


def use_custodian(original_workflow):
    """
    Replaces all tasks with "RunVasp" (e.g. RunVaspDirect) to be
    RunVaspCustodian. Thus, this powerup adds error correction into VASP
    runs if not originally present.

    Args:
        original_workflow (Workflow)
    """
    wf_dict = original_workflow.to_dict()
    vasp_fws_and_tasks = get_runvasp_fws_and_tasks(original_workflow)

    for idx_fw, idx_t in vasp_fws_and_tasks:
        vasp_cmd = wf_dict["fws"][idx_fw]["spec"]["_tasks"][idx_t]["vasp_cmd"]
        wf_dict["fws"][idx_fw]["spec"]["_tasks"][idx_t] = \
            RunVaspCustodian(vasp_cmd=vasp_cmd).to_dict()

    return Workflow.from_dict(wf_dict)


def make_fake_workflow(original_workflow):
    """
    Replaces all tasks with "RunVasp" (e.g. RunVaspDirect) to be
    RunVaspFake. Thus, we do not actually run VASP but copy
    pre-determined inputs and outputs.

    Args:
        original_workflow (Workflow)
    """
    wf_dict = original_workflow.to_dict()
    for idx_fw, fw in enumerate(original_workflow.fws):
        for job_type in fake_dirs.keys():
            if job_type in fw.name:
                for idx_t, t in enumerate(fw.tasks):
                    if is_runvasp_task(t):
                        wf_dict["fws"][idx_fw]["spec"]["_tasks"][
                            idx_t] = RunVaspFake(
                            fake_dir=fake_dirs[job_type]).to_dict()

    return Workflow.from_dict(wf_dict)


def decorate_write_name(original_wf, use_slug=True):
    """
    Every FireWork begins by writing an empty file with the name
    "FW--<fw.name>". This makes it easy to figure out what jobs are in what
    launcher directories, e.g. "ls -l launch*/FW--*" from within a "block" dir.

    Args:
        original_wf (Workflow)
        use_slug (bool): whether to replace whitespace-type chars with a slug
    """
    for fw in original_wf.fws:
        fname = "FW--{}".format(fw.name)
        if use_slug:
            fname = get_slug(fname)
        fw.spec["_tasks"].insert(0, FileWriteTask(
            files_to_write=[{"filename": fname, "contents": ""}]).to_dict())
    return original_wf


def add_trackers(original_wf):
    """
    Every FireWork that runs VASP also tracks the OUTCAR and OSZICAR using FWS Trackers.

    Args:
        original_wf (Workflow)

    """
    tracker1 = Tracker('OUTCAR', nlines=25)
    tracker2 = Tracker('OSZICAR', nlines=25)

    wf_dict = original_wf.to_dict()

    for idx_fw, idx_t in get_runvasp_fws_and_tasks(original_wf):
        if "_trackers" in wf_dict["fws"][idx_fw]["spec"]:
            wf_dict["fws"][idx_fw]["spec"]["_trackers"].extend([tracker1, tracker2])
        else:
            wf_dict["fws"][idx_fw]["spec"]["_trackers"] = [tracker1, tracker2]

    return Workflow.from_dict(wf_dict)