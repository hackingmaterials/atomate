# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import numpy as np

from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet
from pymatgen.core import Structure

from atomate.vasp.config import SMALLGAP_KPOINT_MULTIPLY, STABILITY_CHECK, VASP_CMD, \
    DB_FILE, ADD_WF_METADATA
from atomate.vasp.powerups import add_small_gap_multiply, add_stability_check, add_modify_incar, \
    add_wf_metadata, add_common_powerups
from atomate.vasp.workflows.base.core import get_wf
from atomate.vasp.workflows.base.elastic import get_wf_elastic_constant
from atomate.vasp.workflows.base.raman import get_wf_raman_spectra
from atomate.vasp.workflows.base.gibbs import get_wf_gibbs_free_energy
from atomate.vasp.workflows.base.bulk_modulus import get_wf_bulk_modulus
from atomate.vasp.workflows.base.thermal_expansion import get_wf_thermal_expansion
from atomate.vasp.workflows.base.neb import get_wf_neb_from_endpoints, \
    get_wf_neb_from_structure, get_wf_neb_from_images


__author__ = 'Anubhav Jain, Kiran Mathew'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov'


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

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

    return wf


def wf_bandstructure_plus_hse(structure, gap_only=True, c=None):
    c = c or {}
    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    db_file = c.get("DB_FILE", DB_FILE)

    if gap_only:
        wf_src_name = "bandstructure_hsegap.yaml"
    else:
        wf_src_name = "bandstructure_hse.yaml"

    wf = get_wf(structure, wf_src_name,
                vis=MPRelaxSet(structure, force_gamma=True),
                common_params={"vasp_cmd": vasp_cmd,
                               "db_file": db_file})

    wf = add_common_powerups(wf, c)

    if c.get("SMALLGAP_KPOINT_MULTIPLY", SMALLGAP_KPOINT_MULTIPLY):
        wf = add_small_gap_multiply(wf, 0.5, 5, "static")
        wf = add_small_gap_multiply(wf, 0.5, 5, "nscf")

    if c.get("STABILITY_CHECK", STABILITY_CHECK):
        wf = add_stability_check(wf, fw_name_constraint="structure optimization")

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

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

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

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

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

    return wf


def wf_structure_optimization(structure, c=None):

    c = c or {}
    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    db_file = c.get("DB_FILE", DB_FILE)
    user_incar_settings = c.get("USER_INCAR_SETTINGS")

    wf = get_wf(structure, "optimize_only.yaml",
                vis=MPRelaxSet(structure, force_gamma=True,
                               user_incar_settings=user_incar_settings),
                common_params={"vasp_cmd": vasp_cmd,
                               "db_file": db_file})

    wf = add_common_powerups(wf, c)

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

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

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

    return wf


def wf_piezoelectric_constant(structure, c=None):

    c = c or {}
    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    db_file = c.get("DB_FILE", DB_FILE)

    wf = get_wf(structure, "piezoelectric_constant.yaml",
                common_params={"vasp_cmd": vasp_cmd,
                               "db_file": db_file})

    wf = add_common_powerups(wf, c)

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

    return wf


def wf_elastic_constant(structure, c=None):

    c = c or {}
    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    db_file = c.get("DB_FILE", DB_FILE)
    user_kpoints_settings = c.get("user_kpoints_settings", {"grid_density": 7000})
    norm_deformations = c.get("norm_deformations", [-0.01, -0.005, 0.005, 0.01])
    shear_deformations = c.get("shear_deformations", [-0.06, -0.03, 0.03, 0.06])
    optimize_structure = c.get("optimize_structure", True)

    wf = get_wf_elastic_constant(structure, vasp_cmd=vasp_cmd,
                                 norm_deformations=norm_deformations,
                                 shear_deformations=shear_deformations,
                                 db_file=db_file, user_kpoints_settings=user_kpoints_settings,
                                 optimize_structure=optimize_structure)
    mip = {"incar_update":{"ENCUT": 700, "EDIFF": 1e-6, "LAECHG":False}}
    wf = add_modify_incar(wf, modify_incar_params=mip)

    wf = add_common_powerups(wf, c)

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

    return wf


def wf_raman_spectra(structure, c=None):
    """
    Raman spectra workflow from the given structure and config dict

    Args:
        structure (Structure): input structure
        c (dict): workflow config dict

    Returns:
        Workflow
    """

    c = c or {}
    wf = get_wf_raman_spectra(structure, modes=c.get("modes", None), step_size=c.get("step_size", 0.005),
                              vasp_cmd=c.get("vasp_cmd", "vasp"), db_file=c.get("db_file", None))

    wf = add_modify_incar(wf, modify_incar_params={"incar_update": {"ENCUT": 600, "EDIFF": 1e-6}},
                          fw_name_constraint="static dielectric")

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

    return wf


def wf_gibbs_free_energy(structure, c=None):
    """
    Gibbs free energy workflow from the given structure and config dict.

    Args:
        structure (Structure): input structure
        c (dict): workflow config dict

    Returns:
        Workflow
    """
    c = c or {}

    vasp_cmd = c.get("vasp_cmd", VASP_CMD)
    db_file = c.get("db_file", DB_FILE)
    user_kpoints_settings = c.get("user_kpoints_settings", {"grid_density": 7000})

    # 21 deformed structures
    deformations = c.get("deformations", [(np.identity(3)*(1+x)).tolist()
                                          for x in np.linspace(-0.1, 0.1, 21)])

    eos = c.get("eos", "vinet")
    qha_type = c.get("qha_type", "debye_model")

    # min and max temp in K, default setting is to compute the properties at 300K only
    t_min = c.get("t_min", 300.0)
    t_max = c.get("t_max", 300.0)
    t_step = c.get("t_step", 100.0)

    pressure = c.get("pressure", 0.0)
    poisson = c.get("poisson", 0.25)

    wf = get_wf_gibbs_free_energy(structure, user_kpoints_settings=user_kpoints_settings,
                                  deformations=deformations, vasp_cmd=vasp_cmd, db_file=db_file,
                                  eos=eos, qha_type=qha_type, pressure=pressure, poisson=poisson,
                                  t_min=t_min, t_max=t_max, t_step=t_step)

    wf = add_modify_incar(wf, modify_incar_params={"incar_update": {"ENCUT": 600, "EDIFF": 1e-6}})

    wf = add_common_powerups(wf, c)

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

    return wf


def wf_bulk_modulus(structure, c=None):
    """
    Bulk modulus workflow from the given structure and config dict.

    Args:
        structure (Structure): input structure
        c (dict): workflow config dict

    Returns:
        Workflow
    """
    c = c or {}
    eos = c.get("eos", "vinet")
    vasp_cmd = c.get("vasp_cmd", VASP_CMD)
    db_file = c.get("db_file", DB_FILE)
    user_kpoints_settings = c.get("user_kpoints_settings", {"grid_density": 7000})
    deformations = c.get("deformations", [(np.identity(3)*(1+x)).tolist()
                                          for x in np.linspace(-0.1, 0.1, 10)])

    wf = get_wf_bulk_modulus(structure, eos=eos, user_kpoints_settings=user_kpoints_settings,
                             deformations=deformations, vasp_cmd=vasp_cmd, db_file=db_file)

    wf = add_modify_incar(wf, modify_incar_params={"incar_update": {"ENCUT": 600, "EDIFF": 1e-6}})

    wf = add_common_powerups(wf, c)

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

    return wf


def wf_thermal_expansion(structure, c=None):
    """
    Thermal expansion coefficient workflow from the given structure and config dict.

    Args:
        structure (Structure): input structure
        c (dict): workflow config dict

    Returns:
        Workflow
    """
    c = c or {}
    eos = c.get("eos", "vinet")
    vasp_cmd = c.get("vasp_cmd", VASP_CMD)
    db_file = c.get("db_file", DB_FILE)
    user_kpoints_settings = c.get("user_kpoints_settings", {"grid_density": 7000})
    deformations = c.get("deformations", [(np.identity(3)*(1+x)).tolist()
                                          for x in np.linspace(-0.1, 0.1, 10)])
    pressure = c.get("pressure", 0.0)

    wf = get_wf_thermal_expansion(structure, user_kpoints_settings=user_kpoints_settings,
                                  deformations=deformations, vasp_cmd=vasp_cmd, db_file=db_file,
                                  eos=eos, pressure=pressure)

    wf = add_modify_incar(wf, modify_incar_params={"incar_update": {"ENCUT": 600, "EDIFF": 1e-6}})

    wf = add_common_powerups(wf, c)

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

    return wf


def wf_nudged_elastic_band(structures, c=None, parent=None):
    """
    Nudged elastic band workflow from the given structures and config dict.

    'is_optimized' = config, otherwise False
    'neb_round' = config, otherwise 1

    Args:
        structures (Structure / [Structure]): input structures
        c (dict): workflow config dict, basic format:
                      {"fireworks": [], "is_optimized": False, "common_params": {}}
        parent (Structure): parent structure when not given.
    Returns:
        Workflow
    """
    def get_incar(mode):
        """Get user_incar_settings for fireworks."""

        uis_ini, uis_ep, uis_neb = {}, {}, [{}] * neb_round

        assert mode in [1, 2, 3, 4, 5], "Unknown mode!"
        if mode == 1:
            uis_ini = c["fireworks"][0].get("user_incar_settings", {})
            uis_ep = c["fireworks"][1].get("user_incar_settings", {})
            for i in range(neb_round):
                try:
                    uis_neb[i] = c["fireworks"][3 + i]["user_incar_settings"]
                except:
                    uis_neb[i] = {}
        elif mode in [2, 3]:
            uis_ep = c["fireworks"][0].get("user_incar_settings", {})
            for i in range(neb_round):
                try:
                    uis_neb[i] = c["fireworks"][2 + i]["user_incar_settings"]
                except:
                    uis_neb[i] = {}
        else:
            for i in range(neb_round):
                try:
                    uis_neb[i] = c["fireworks"][i]["user_incar_settings"]
                except:
                    uis_neb[i] = {}
        return uis_ini, uis_ep, uis_neb

    # Structure --> [Structure]
    if isinstance(structures, Structure):
        structures = [structures]
    if not(isinstance(structures, list) and len(structures) > 0):
        raise ValueError("structures must be a list of Structure!")

    # config initialization
    c = c or {}
    spec = c.get("common_params", {})  # default
    is_optimized = c.get("is_optimized", False)  # default
    wf_name = spec.get("wf_name")  # default
    path_sites = spec.get("path_sites", [])  # default

    # Default
    if "fireworks" not in c:
        neb_round = 1
        if len(structures) == 1:
            mode = 2 if is_optimized else 1
        elif len(structures) == 2:
            mode = 4 if is_optimized else 3
        else:
            mode = 5
    # construct wf using config file
    else:
        fw_list = [f['fw'] for f in c.get("fireworks")]
        neb_round = fw_list.count("atomate.vasp.fireworks.core.NEBFW")
        assert neb_round > 0, "No NEB fireworks in config file!"
        # Get mode number
        if len(structures) == 1:
            if not is_optimized and len(fw_list) - neb_round == 3:
                mode = 1
            elif is_optimized and len(fw_list) - neb_round == 2:
                mode = 2
            else:
                raise ValueError("structure conflict with config file settings!")
        elif len(structures) == 2:
            assert isinstance(parent, Structure), "Parent structure is not provided!"
            if not is_optimized and len(fw_list) - neb_round == 2:
                mode = 3
            elif is_optimized and len(fw_list) - neb_round == 0:
                mode = 4
            else:
                raise ValueError("structure conflict with config file settings!")
        else:  # len(structures) >= 3
            assert isinstance(parent, Structure), "Parent structure is not provided!"
            mode = 5

    uis_ini, uis_ep, uis_neb = get_incar(mode)

    user_incar_settings = {"parent": uis_ini, "endpoints": uis_ep, "NEB": uis_neb}

    # Assign workflow using mode
    if mode in [1, 2]:
        wf = get_wf_neb_from_structure(structure=structures[0], path_sites=path_sites,
                                       is_optimized=is_optimized, wf_name=wf_name,
                                       spec=spec, user_incar_settings=user_incar_settings)
    elif mode in [3, 4]:
        wf = get_wf_neb_from_endpoints(endpoints=structures, is_optimized=is_optimized,
                                       wf_name=wf_name, additional_spec=spec,
                                       user_incar_settings=user_incar_settings)
    else:  # mode == 5
        wf = get_wf_neb_from_images(parent, user_incar_settings, images=structures,
                                    wf_name=wf_name, spec=spec)
    return wf