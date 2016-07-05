from __future__ import absolute_import

from matmethods.vasp.vasp_config import ADD_NAMEFILE, SCRATCH_DIR, \
    SMALLGAP_KPOINT_MULTIPLY, ADD_MODIFY_INCAR, STABILITY_CHECK, VASP_CMD, \
    DB_FILE
from matmethods.vasp.workflows.base.core import get_wf
from matmethods.vasp.vasp_powerups import add_namefile, \
    add_small_gap_multiply, use_scratch_dir, add_stability_check, \
    add_modify_incar
from matmethods.vasp.workflows.base.elastic import get_wf_elastic_constant
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


# TODO: config dict stuff could be clearer

def wf_bandstructure(structure, c=None):
    c = c or {}
    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    db_file = c.get("DB_FILE", DB_FILE)

    wf = get_wf(structure, "bandstructure.yaml",
                vis=MPRelaxSet(structure, force_gamma=True),
                common_params={"vasp_cmd": vasp_cmd,
                               "db_file": db_file})

    wf = add_common_powerups(wf, c)

    if c.get("SMALLGAP_KPOINT_MULTIPLY", SMALLGAP_KPOINT_MULTIPLY):
        wf = add_small_gap_multiply(wf, 0.5, 5, "static")
        wf = add_small_gap_multiply(wf, 0.5, 5, "nscf")

    if c.get("STABILITY_CHECK", STABILITY_CHECK):
        wf = add_stability_check(wf, fw_name_constraint="structure optimization")

    return wf


def wf_bandstructure_plus_hse(structure, c=None):
    c = c or {}
    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    db_file = c.get("DB_FILE", DB_FILE)

    wf = get_wf(structure, "bandstructure_hsegap.yaml",
                vis=MPRelaxSet(structure, force_gamma=True),
                common_params={"vasp_cmd": vasp_cmd,
                               "db_file": db_file})

    wf = add_common_powerups(wf, c)

    if c.get("SMALLGAP_KPOINT_MULTIPLY", SMALLGAP_KPOINT_MULTIPLY):
        wf = add_small_gap_multiply(wf, 0.5, 5, "static")
        wf = add_small_gap_multiply(wf, 0.5, 5, "nscf")

    if c.get("STABILITY_CHECK", STABILITY_CHECK):
        wf = add_stability_check(wf, fw_name_constraint="structure optimization")

    return wf


def wf_bandstructure_plus_boltztrap(structure, c=None):
    c = c or {}
    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    db_file = c.get("DB_FILE", DB_FILE)

    params = []
    for x in range(4):
        params.append({"vasp_cmd": vasp_cmd, "db_file": db_file})
    params.append({"db_file": db_file})

    wf = get_wf(structure, "bandstructure_boltztrap.yaml",
                vis=MPRelaxSet(structure, force_gamma=True),
                params=params)

    wf = add_common_powerups(wf, c)

    if c.get("SMALLGAP_KPOINT_MULTIPLY", SMALLGAP_KPOINT_MULTIPLY):
        wf = add_small_gap_multiply(wf, 0.5, 5, "static")
        wf = add_small_gap_multiply(wf, 0.5, 5, "nscf")

    if c.get("STABILITY_CHECK", STABILITY_CHECK):
        wf = add_stability_check(wf, fw_name_constraint="structure optimization")

    return wf


def wf_static(structure, c=None):
    c = c or {}
    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    db_file = c.get("DB_FILE", DB_FILE)

    wf = get_wf(structure, "static_only.yaml",
                vis=MPStaticSet(structure),
                common_params={"vasp_cmd": vasp_cmd,
                               "db_file": db_file})

    wf = add_common_powerups(wf, c)
    return wf


def wf_structure_optimization(structure, c=None):

    c = c or {}
    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    db_file = c.get("DB_FILE", DB_FILE)

    wf = get_wf(structure, "optimize_only.yaml",
                vis=MPRelaxSet(structure, force_gamma=True),
                common_params={"vasp_cmd": vasp_cmd,
                               "db_file": db_file})

    wf = add_common_powerups(wf, c)
    return wf


def wf_dielectric_constant(structure, c=None):

    c = c or {}
    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    db_file = c.get("DB_FILE", DB_FILE)

    wf = get_wf(structure, "dielectric_constant.yaml",
                vis=MPRelaxSet(structure, force_gamma=True),
                common_params={"vasp_cmd": vasp_cmd,
                               "db_file": db_file})

    wf = add_common_powerups(wf, c)
    return wf


def wf_piezoelectric_constant(structure, c=None):

    c = c or {}
    wf = wf_dielectric_constant(structure, c)

    wf = add_modify_incar(wf, modify_incar_params={"incar_update": {"ENCUT": 1000,
                                                                    "ADDGRID": True,
                                                                    "LREAL": False,
                                                                    "EDIFF": 1e-7}
                                                   },
                          fw_name_constraint="static dielectric")
    for fw in wf.fws:
        fw.name = fw.name.replace("dielectric", "piezoelectric")

    return wf


def wf_elastic_constant(structure, c=None):

    c = c or {}
    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    db_file = c.get("DB_FILE", DB_FILE)

    wf = get_wf_elastic_constant(structure, vasp_cmd=vasp_cmd,
                                 db_file=db_file)

    wf = add_modify_incar(wf, modify_incar_params={"incar_update":
                                                   {"ENCUT": 700, "EDIFF": 1e-6}})

    wf = add_common_powerups(wf, c)
    return wf


def add_common_powerups(wf, c):

    if c.get("ADD_NAMEFILE", ADD_NAMEFILE):
        wf = add_namefile(wf)

    if c.get("SCRATCH_DIR", SCRATCH_DIR):
        wf = use_scratch_dir(wf, c.get("SCRATCH_DIR", SCRATCH_DIR))

    if c.get("ADD_MODIFY_INCAR", ADD_MODIFY_INCAR):
        wf = add_modify_incar(wf)

    return wf
