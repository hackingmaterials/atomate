# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

from pymatgen.io.vasp.sets import MVLScanRelaxSet

from atomate.vasp.config import VASP_CMD, DB_FILE, ADD_WF_METADATA, HALF_KPOINTS_FIRST_RELAX, REMOVE_WAVECAR
from atomate.vasp.powerups import use_custodian, add_wf_metadata, add_common_powerups, clean_up_files
from atomate.vasp.workflows.base.core import get_wf

__author__ = 'Shyam Dwaraknath, Anubhav Jain'
__email__ = 'shyamd@lbl.gov, ajain@lbl.gov'


def wf_scan_opt(structure, c=None):

    c = c or {}
    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    db_file = c.get("DB_FILE", DB_FILE)
    user_incar_settings = c.get("USER_INCAR_SETTINGS", {})
    half_kpts = c.get("HALF_KPOINTS_FIRST_RELAX", HALF_KPOINTS_FIRST_RELAX)
    ediffg = user_incar_settings.get("EDIFFG", -0.05)

    wf = get_wf(
        structure,
        "optimize_only.yaml",
        vis=MVLScanRelaxSet(
            structure, user_incar_settings=user_incar_settings),
        common_params={
            "vasp_cmd": vasp_cmd,
            "db_file": db_file
        })

    wf = use_custodian(
        wf,
        custodian_params={
            "ediffg": ediffg,
            "max_force_threshold": 0,
            "half_kpts_first_relax": half_kpts,
            "job_type": "metagga_opt_run",
            "db_file": db_file,
            "vasp_cmd": vasp_cmd
        })
    wf = add_common_powerups(wf, c)

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

    if c.get("REMOVE_WAVECAR", REMOVE_WAVECAR):
        wf = clean_up_files(wf)

    return wf
