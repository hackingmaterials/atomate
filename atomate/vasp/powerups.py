# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

from fireworks import Workflow, FileWriteTask
from fireworks.core.firework import Tracker
from fireworks.utilities.fw_utilities import get_slug

from pymatgen import Structure

from atomate.utils.utils import get_meta_from_structure, get_fws_and_tasks
from atomate.vasp.firetasks.glue_tasks import CheckStability, CheckBandgap
from atomate.vasp.firetasks.run_calc import RunVaspCustodian, RunVaspFake, RunVaspDirect, RunNoVasp
from atomate.vasp.firetasks.neb_tasks import RunNEBVaspFake
from atomate.vasp.firetasks.write_inputs import ModifyIncar
from atomate.vasp.firetasks.parse_outputs import JsonToDb
from atomate.vasp.config import ADD_NAMEFILE, SCRATCH_DIR, ADD_MODIFY_INCAR, GAMMA_VASP_CMD

__author__ = 'Anubhav Jain, Kiran Mathew'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov'


def add_priority(original_wf, root_priority, child_priority=None):
    """
    Adds priority to a workflow

    Args:
        original_wf (Workflow): original WF
        root_priority (int): priority of first (root) job(s)
        child_priority(int): priority of all child jobs. Defaults to root_priority

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


def remove_custodian(original_wf, fw_name_constraint=None):
    """
    Replaces all tasks with "RunVasp*" (e.g. RunVaspCustodian) to be RunVaspDirect.

    Args:
        original_wf (Workflow): original workflow
        fw_name_constraint (str): Only apply changes to FWs where fw_name contains this substring.

    Returns:
       Workflow
    """
    vasp_fws_and_tasks = get_fws_and_tasks(original_wf, fw_name_constraint=fw_name_constraint,
                                           task_name_constraint="RunVasp")
    for idx_fw, idx_t in vasp_fws_and_tasks:
        vasp_cmd = original_wf.fws[idx_fw].tasks[idx_t]["vasp_cmd"]
        original_wf.fws[idx_fw].tasks[idx_t] = RunVaspDirect(vasp_cmd=vasp_cmd)
    return original_wf


def use_custodian(original_wf, fw_name_constraint=None, custodian_params=None):
    """
    Replaces all tasks with "RunVasp*" (e.g. RunVaspDirect) to be RunVaspCustodian. Thus, this
    powerup adds error correction into VASP runs if not originally present and/or modifies
    the correction behavior.

    Args:
        original_wf (Workflow): original workflow
        fw_name_constraint (str): Only apply changes to FWs where fw_name contains this substring.
            For example, use custodian only for certain runs, or set job_type to
            "double_relaxation_run" only for structure optimization run, or set different
            handler_group for different runs.
        custodian_params (dict): A dict of parameters for RunVaspCustodian. e.g., use it to set
            a "scratch_dir" or "handler_group".

    Returns:
       Workflow
    """
    custodian_params = custodian_params if custodian_params else {}
    vasp_fws_and_tasks = get_fws_and_tasks(original_wf, fw_name_constraint=fw_name_constraint,
                                           task_name_constraint="RunVasp")
    for idx_fw, idx_t in vasp_fws_and_tasks:
        if "vasp_cmd" not in custodian_params:
            custodian_params["vasp_cmd"] = original_wf.fws[idx_fw].tasks[idx_t]["vasp_cmd"]
        original_wf.fws[idx_fw].tasks[idx_t] = \
            RunVaspCustodian(**custodian_params)
    return original_wf


def use_no_vasp(original_wf, ref_dirs):
    """
    Instead of running VASP, does nothing and pass task documents from task.json files in ref_dirs
    to task database.

    Args:
        original_wf (Workflow)
        ref_dirs(dict): key=firework name, value=path to the reference vasp calculation directory

    Returns:
        Workflow
    """
    for idx_fw, fw in enumerate(original_wf.fws):
        for job_type in ref_dirs.keys():
            if job_type in fw.name:
                for idx_t, t in enumerate(fw.tasks):
                    if "RunVasp" in str(t):
                        original_wf.fws[idx_fw].tasks[idx_t] = RunNoVasp(ref_dir=ref_dirs[job_type])
                    if "VaspToDb" in str(t):
                        original_wf.fws[idx_fw].tasks[idx_t] = JsonToDb(db_file=t.get("db_file", None),
                                                                        calc_dir=ref_dirs[job_type])
    return original_wf


def use_fake_vasp(original_wf, ref_dirs, params_to_check=None):
    """
    Replaces all tasks with "RunVasp" (e.g. RunVaspDirect) to be RunVaspFake. Thus, we do not
    actually run VASP but copy pre-determined inputs and outputs.

    Args:
        original_wf (Workflow)
        ref_dirs (dict): key=firework name, value=path to the reference vasp calculation directory
        params_to_check (list): optional list of incar parameters to check.

    Returns:
        Workflow
    """
    if not params_to_check:
        params_to_check = ["ISPIN", "ENCUT", "ISMEAR", "SIGMA", "IBRION", "LORBIT",
                           "NBANDS", "LMAXMIX"]
    for idx_fw, fw in enumerate(original_wf.fws):
        for job_type in ref_dirs.keys():
            if job_type in fw.name:
                for idx_t, t in enumerate(fw.tasks):
                    if "RunVasp" in str(t):
                        original_wf.fws[idx_fw].tasks[idx_t] = \
                            RunVaspFake(ref_dir=ref_dirs[job_type],
                                        params_to_check=params_to_check)
                    if "RunVaspCustodian" in str(t) and t.get("job_type") == "neb":
                        original_wf.fws[idx_fw].tasks[idx_t] = \
                            RunNEBVaspFake(ref_dir=ref_dirs[job_type],
                                           params_to_check=params_to_check)

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
        original_wf.fws[idx].tasks.insert(0, FileWriteTask(
            files_to_write=[{"filename": fname, "contents": ""}]))
    return original_wf


def add_trackers(original_wf, tracked_files=None, nlines=25):
    """
    Every FireWork that runs VASP also tracks the OUTCAR, OSZICAR, etc using FWS Trackers.

    Args:
        original_wf (Workflow)
        tracked_files (list) : list of files to be tracked
        nlines (int): number of lines at the end of files to be tracked

    Returns:
       Workflow
    """
    if tracked_files is None:
        tracked_files = ["OUTCAR", "OSZICAR"]
    trackers = [Tracker(f, nlines=nlines, allow_zipped=True) for f in tracked_files]
    for idx_fw, idx_t in get_fws_and_tasks(original_wf, task_name_constraint="RunVasp"):
        if "_trackers" in original_wf.fws[idx_fw].spec:
            original_wf.fws[idx_fw].spec["_trackers"].extend(trackers)
        else:
            original_wf.fws[idx_fw].spec["_trackers"] = trackers
    return original_wf


def add_modify_incar(original_wf, modify_incar_params=None, fw_name_constraint=None):
    """
    Every FireWork that runs VASP has a ModifyIncar task just beforehand. For example, allows
    you to modify the INCAR based on the Worker using env_chk or using hard-coded changes.

    Args:
        original_wf (Workflow)
        modify_incar_params (dict) - dict of parameters for ModifyIncar.
        fw_name_constraint (str) - Only apply changes to FWs where fw_name contains this substring.

    Returns:
       Workflow
    """
    modify_incar_params = modify_incar_params or {"incar_update": ">>incar_update<<"}
    for idx_fw, idx_t in get_fws_and_tasks(original_wf, fw_name_constraint=fw_name_constraint,
                                           task_name_constraint="RunVasp"):
        original_wf.fws[idx_fw].tasks.insert(idx_t, ModifyIncar(**modify_incar_params))
    return original_wf


def modify_to_soc(original_wf, nbands, structure=None, modify_incar_params=None, fw_name_constraint=None):
    """
    Takes a regular workflow and transforms its VASP fireworkers that are specified with
    fw_name_constraints to non-collinear calculations taking spin orbit coupling into account.

    Args:
        original_wf (Workflow)
        nbands (int): number of bands selected by the user (for now)
        structure (Structure)
        modify_incar_params ({}): a dictionary containing the setting for modyfining the INCAR (e.g. {"ICHARG": 11})
        fw_name_constraint (string): name of the fireworks to be modified (all if None is passed)

    Returns:
        Workflow: modified with SOC
    """

    if structure is None:
        try:
            sid = get_fws_and_tasks(original_wf, fw_name_constraint="structure optimization",
                                    task_name_constraint="RunVasp")[0][0]
            structure = Structure.from_dict(original_wf.fws[sid].tasks[1]["vasp_input_set"]["structure"])
        except:
            raise ValueError("modify_to_soc powerup requires the structure in vasp_input_set")

    magmom = ""
    for _ in structure:
        magmom += "0 0 0.6 "
    #TODO: add saxis as an input parameter with default being (0 0 1)
    modify_incar_params = modify_incar_params or {"incar_update": {"LSORBIT": "T", "NBANDS":
        nbands, "MAGMOM": magmom, "ISPIN": 1, "LMAXMIX": 4, "ISYM": 0}}

    for idx_fw, idx_t in get_fws_and_tasks(original_wf, fw_name_constraint=fw_name_constraint,
                                           task_name_constraint="RunVasp"):
        if "structure" not in original_wf.fws[idx_fw].name and "static" not in original_wf.fws[idx_fw].name:
            original_wf.fws[idx_fw].tasks[idx_t]["vasp_cmd"] = ">>vasp_ncl<<"
            original_wf.fws[idx_fw].tasks.insert(idx_t, ModifyIncar(**modify_incar_params))

        original_wf.fws[idx_fw].name += " soc"

    for idx_fw, idx_t in get_fws_and_tasks(original_wf, fw_name_constraint=fw_name_constraint,
                                           task_name_constraint="RunBoltztrap"):
        original_wf.fws[idx_fw].name += " soc"

    return original_wf


def set_fworker(original_wf, fworker_name, fw_name_constraint=None, task_name_constraint=None):
    """
    set _fworker spec of Fireworker(s) of a Workflow. It can be used to specify a queue;
    e.g. run large-memory jobs on a separate queue.

    Args:
        original_wf (Workflow):
        fworker_name (str): user-defined tag to be added under fw.spec._fworker
            e.g. "large memory", "big", etc
        fw_name_constraint (str): name of the Fireworks to be tagged (all if None is passed)
        task_name_constraint (str): name of the Firetasks to be tagged (e.g. None or 'RunVasp')

    Returns:
        Workflow: modified workflow with specified Fireworkers tagged
    """
    for idx_fw, idx_t in get_fws_and_tasks(original_wf, fw_name_constraint=fw_name_constraint,
                                           task_name_constraint=task_name_constraint):
        original_wf.fws[idx_fw].spec["_fworker"] = fworker_name
    return original_wf


def add_wf_metadata(original_wf, structure):
    """
    Adds structure metadata to a workflow

    Args:
        original_wf: (Workflow)
        structure: (Structure) the structure being run by this workflow

    Returns:
        Workflow
    """
    original_wf.metadata["structure"] = structure.as_dict()
    original_wf.metadata.update(get_meta_from_structure(structure))
    return original_wf


def add_stability_check(original_wf, check_stability_params=None, fw_name_constraint=None):
    """
    Every FireWork that enters into the Db has a CheckStability task afterward. This
    allows defusing jobs that are not stable. In practice, you might want
    to set the fw_name_constraint so that the stability is only checked at the
    beginning of the workflow

    Args:
        original_wf (Workflow)
        check_stability_params (dict): a **kwargs** style dict of params
        fw_name_constraint (str) - Only apply changes to FWs where fw_name contains this substring.

    Returns:
       Workflow
    """
    check_stability_params = check_stability_params or {}
    for idx_fw, idx_t in get_fws_and_tasks(original_wf, fw_name_constraint=fw_name_constraint,
                                           task_name_constraint="VaspToDb"):
        original_wf.fws[idx_fw].tasks.append(CheckStability(**check_stability_params))
    return original_wf


def add_bandgap_check(original_wf, check_bandgap_params=None, fw_name_constraint=None):
    """
    Every FireWork that enters into the Db has a band gap check afterwards, e.g. min_gap and max_gap

    Args:
        original_wf (Workflow)
        check_bandgap_params (dict): a **kwargs** style dict of params, e.g. min_gap or max_gap
        fw_name_constraint (str) - Only apply changes to FWs where fw_name contains this substring.

    Returns:
       Workflow
    """
    check_bandgap_params = check_bandgap_params or {}
    for idx_fw, idx_t in get_fws_and_tasks(original_wf, fw_name_constraint=fw_name_constraint,
                                           task_name_constraint="VaspToDb"):
        original_wf.fws[idx_fw].tasks.append(CheckBandgap(**check_bandgap_params))
    return original_wf


def add_modify_incar_envchk(original_wf, fw_name_constraint=None):
    """
    If you set the "incar_update" parameter in the Worker env, the INCAR will update this
    parameter for all matching VASP runs

    Args:
        original_wf (Workflow)
        fw_name_constraint (str) - Only apply changes to FWs where fw_name contains this substring.

    Returns:
       Workflow
    """
    return add_modify_incar(original_wf, {"incar_update": ">>incar_update<<"},
                            fw_name_constraint=fw_name_constraint)


def add_small_gap_multiply(original_wf, gap_cutoff, density_multiplier, fw_name_constraint=None):
    """
    In all FWs with specified name constraints, add a 'small_gap_multiply' parameter that
    multiplies the k-mesh density of compounds with gap < gap_cutoff by density multiplier.
    Useful for increasing the k-point mesh for metallic or small gap systems.
    Note that this powerup only works on FireWorks with the appropriate WriteVasp* tasks that
    accept the small_gap_multiply argument...

    Args:
        original_wf (Workflow)
        gap_cutoff (float): Only multiply k-points for materials with gap < gap_cutoff (eV)
        density_multiplier (float): Multiply k-point density by this amount
        fw_name_constraint (str): Only apply changes to FWs where fw_name contains this substring.

    Returns:
       Workflow
    """
    for idx_fw, idx_t in get_fws_and_tasks(original_wf, fw_name_constraint=fw_name_constraint,
                                           task_name_constraint="WriteVasp"):
        original_wf.fws[idx_fw].tasks[idx_t]["small_gap_multiply"] = [gap_cutoff, density_multiplier]
    return original_wf


def use_scratch_dir(original_wf, scratch_dir):
    """
    For all RunVaspCustodian tasks, add the desired scratch dir.

    Args:
        original_wf (Workflow)
        scratch_dir (path): Path to the scratch dir to use. Supports env_chk

    Returns:
       Workflow
    """
    for idx_fw, idx_t in get_fws_and_tasks(original_wf, task_name_constraint="RunVaspCustodian"):
        original_wf.fws[idx_fw].tasks[idx_t]["scratch_dir"] = scratch_dir
    return original_wf


def add_additional_fields_to_taskdocs(original_wf, update_dict=None):
    """
    For all VaspToDbTasks in a given workflow, add information  to "additional_fields" to be
    placed in the task doc.
    
    Args:
        original_wf (Workflow)
        update_dict (Dict): dictionary to add additional_fields

    Returns:
       Workflow
    """
    for idx_fw, idx_t in get_fws_and_tasks(original_wf, task_name_constraint="VaspToDb"):
        original_wf.fws[idx_fw].tasks[idx_t]["additional_fields"].update(update_dict)
    return original_wf


def add_tags(original_wf, tags_list):
    """
    Adds tags to all Fireworks in the Workflow, WF metadata, as well as additional_fields for
    the VaspDrone to track them later (e.g. all fireworks and vasp tasks related to a research project)

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
    for constraint in ["VaspToDb", "BoltztrapToDb"]:
        for idx_fw, idx_t in get_fws_and_tasks(original_wf, task_name_constraint=constraint):
            if "tags" in original_wf.fws[idx_fw].tasks[idx_t]["additional_fields"]:
                original_wf.fws[idx_fw].tasks[idx_t]["additional_fields"][
                    "tags"].extend(tags_list)
            else:
                original_wf.fws[idx_fw].tasks[idx_t]["additional_fields"][
                    "tags"] = tags_list

    return original_wf


def add_common_powerups(wf, c):
    """
    Apply the common powerups such as add_namefile, use_scratch_dir etc. from the given config dict.

    Args:
        wf (Workflow)
        c (dict): Config dict

    Returns:
        Workflow
    """
    if c.get("ADD_NAMEFILE", ADD_NAMEFILE):
        wf = add_namefile(wf)

    if c.get("SCRATCH_DIR", SCRATCH_DIR):
        wf = use_scratch_dir(wf, c.get("SCRATCH_DIR", SCRATCH_DIR))

    if c.get("ADD_MODIFY_INCAR", ADD_MODIFY_INCAR):
        wf = add_modify_incar(wf)

    if c.get("GAMMA_VASP_CMD", GAMMA_VASP_CMD):
        wf = use_gamma_vasp((wf),c.get("GAMMA_VASP_CMD", GAMMA_VASP_CMD))

    return wf


def use_gamma_vasp(original_wf, gamma_vasp_cmd):
    """
    For all RunVaspCustodian tasks, add the desired scratch dir.

    Args:
        original_wf (Workflow)
        gamma_vasp_cmd (str): path to gamma_vasp_cmd. Supports env_chk

    Returns:
       Workflow
    """
    for idx_fw, idx_t in get_fws_and_tasks(original_wf, task_name_constraint="RunVaspCustodian"):
        original_wf.fws[idx_fw].tasks[idx_t]["gamma_vasp_cmd"] = gamma_vasp_cmd
    return original_wf
