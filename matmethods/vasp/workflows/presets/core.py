from __future__ import absolute_import

from matmethods.vasp.workflows.base.core import get_wf

from matmethods.vasp.vasp_powerups import add_namefile, \
    add_small_gap_multiply, use_scratch_dir, add_stability_check, \
    add_modify_incar
from matmethods.vasp.workflows.base.elastic import get_wf_elastic_constant
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

# TODO: clean up some code duplication in config params
# TODO: this needs massive duplication cleanup - simple but just need to do it
def wf_bandstructure(structure, config=None):
    config = config or {}
    wf = get_wf(structure, "bandstructure.yaml",
                vis=MPRelaxSet(structure, force_gamma=True),
                common_params={"vasp_cmd": ">>vasp_cmd<<",
                               "db_file": ">>db_file<<"})

    if config.get("ADD_NAMEFILE", True):
        wf = add_namefile(wf)

    if config.get("SMALLGAP_KPOINT_MULTIPLY", True):
        wf = add_small_gap_multiply(wf, 0.5, 5, "static")
        wf = add_small_gap_multiply(wf, 0.5, 5, "nscf")

    if config.get("USE_SCRATCH_DIR", True):
        wf = use_scratch_dir(wf, ">>scratch_dir<<")

    if config.get("ADD_MODIFY_INCAR", False):
        wf = add_modify_incar(wf)

    if config.get("CHECK_STABILITY", False):
        wf = add_stability_check(wf, fw_name_constraint="structure optimization")
    return wf


def wf_bandstructure_plus_hse(structure, config=None):
    config = config or {}

    wf = get_wf(structure, "bandstructure_hsegap.yaml",
                vis=MPRelaxSet(structure, force_gamma=True),
                common_params={"vasp_cmd": ">>vasp_cmd<<",
                               "db_file": ">>db_file<<"})

    if config.get("ADD_NAMEFILE", True):
        wf = add_namefile(wf)

    if config.get("SMALLGAP_KPOINT_MULTIPLY", True):
        wf = add_small_gap_multiply(wf, 0.5, 5, "static")
        wf = add_small_gap_multiply(wf, 0.5, 5, "nscf")

    if config.get("USE_SCRATCH_DIR", True):
        wf = use_scratch_dir(wf, ">>scratch_dir<<")

    if config.get("ADD_MODIFY_INCAR", False):
        wf = add_modify_incar(wf)

    if config.get("CHECK_STABILITY", False):
        wf = add_stability_check(wf, fw_name_constraint="structure optimization")

    return wf


def wf_bandstructure_plus_boltztrap(structure, config=None):
    config = config or {}

    params = []
    for x in range(4):
        params.append({"vasp_cmd": ">>vasp_cmd<<", "db_file": ">>db_file<<"})
    params.append({})  # no params for Boltztrap FW

    wf = get_wf(structure, "bandstructure_boltztrap.yaml",
                vis=MPRelaxSet(structure, force_gamma=True),
                params=params)

    if config.get("ADD_NAMEFILE", True):
        wf = add_namefile(wf)

    if config.get("SMALLGAP_KPOINT_MULTIPLY", True):
        wf = add_small_gap_multiply(wf, 0.5, 5, "static")
        wf = add_small_gap_multiply(wf, 0.5, 5, "nscf")

    if config.get("USE_SCRATCH_DIR", True):
        wf = use_scratch_dir(wf, ">>scratch_dir<<")

    if config.get("ADD_MODIFY_INCAR", False):
        wf = add_modify_incar(wf)

    if config.get("CHECK_STABILITY", False):
        wf = add_stability_check(wf, fw_name_constraint="structure optimization")

    return wf


def wf_static(structure, config=None):
    config = config or {}

    wf = get_wf(structure, "static_only.yaml",
                vis=MPStaticSet(structure),
                common_params={"vasp_cmd": ">>vasp_cmd<<",
                               "db_file": ">>db_file<<"})

    if config.get("ADD_NAMEFILE", True):
        wf = add_namefile(wf)

    if config.get("ADD_MODIFY_INCAR", False):
        wf = add_modify_incar(wf)

    if config.get("USE_SCRATCH_DIR", True):
        wf = use_scratch_dir(wf, ">>scratch_dir<<")

    return wf


def wf_structure_optimization(structure, config=None):

    config = config or {}
    wf = get_wf(structure, "optimize_only.yaml",
                vis=MPRelaxSet(structure, force_gamma=True),
                common_params={"vasp_cmd": ">>vasp_cmd<<",
                               "db_file": ">>db_file<<"})

    if config.get("ADD_NAMEFILE", True):
        wf = add_namefile(wf)

    if config.get("ADD_MODIFY_INCAR", False):
        wf = add_modify_incar(wf)

    if config.get("USE_SCRATCH_DIR", True):
        wf = use_scratch_dir(wf, ">>scratch_dir<<")

    return wf


def wf_dielectric_constant(structure, config=None):

    config = config or {}

    wf = get_wf(structure, "dielectric_constant.yaml",
                vis=MPRelaxSet(structure, force_gamma=True),
                common_params={"vasp_cmd": ">>vasp_cmd<<",
                               "db_file": ">>db_file<<"})

    if config.get("ADD_NAMEFILE", True):
        wf = add_namefile(wf)

    if config.get("ADD_MODIFY_INCAR", False):
        wf = add_modify_incar(wf)

    if config.get("USE_SCRATCH_DIR", True):
        wf = use_scratch_dir(wf, ">>scratch_dir<<")

    return wf


def wf_piezoelectric_constant(structure, config=None):
    wf = wf_dielectric_constant(structure, config)

    wf = add_modify_incar(wf, modify_incar_params={"incar_update": {"ENCUT": 1000,
                                                                    "ADDGRID": True,
                                                                    "LREAL": False,
                                                                    "EDIFF": 1e-7}
                                                   },
                          fw_name_constraint="static dielectric")
    for fw in wf.fws:
        fw.name = fw.name.replace("dielectric", "piezoelectric")

    return wf


def wf_elastic_constant(structure, config=None):
    wf = get_wf_elastic_constant(structure, vasp_cmd=">>vasp_cmd<<",
                                 db_file=">>db_file<<")

    wf = add_modify_incar(wf, modify_incar_params={"incar_update":
                                                       {"ENCUT": 700, "EDIFF": 1e-6}})
    return wf
                         

