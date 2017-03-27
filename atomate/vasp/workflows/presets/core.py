# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import numpy as np

from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet

from atomate.vasp.config import SMALLGAP_KPOINT_MULTIPLY, STABILITY_CHECK, VASP_CMD, DB_FILE, \
    ADD_WF_METADATA
from atomate.vasp.powerups import add_small_gap_multiply, add_stability_check, add_modify_incar, \
    add_wf_metadata, add_common_powerups
from atomate.vasp.workflows.base.core import get_wf
from atomate.vasp.workflows.base.elastic import get_wf_elastic_constant
from atomate.vasp.workflows.base.raman import get_wf_raman_spectra
from atomate.vasp.workflows.base.gibbs import get_wf_gibbs_free_energy
from atomate.vasp.workflows.base.bulk_modulus import get_wf_bulk_modulus
from atomate.vasp.workflows.base.thermal_expansion import get_wf_thermal_expansion
from atomate.vasp.workflows.base.neb import get_wf_neb_from_endpoints, get_wf_neb_from_structure, \
    get_wf_neb_from_images


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


def wf_dielectric_constant_no_opt(structure, c=None):

    c = c or {}
    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    db_file = c.get("DB_FILE", DB_FILE)

    wf = get_wf(structure, "dielectric_constant_no_opt.yaml",
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
    metadata = c.get("metadata", None)

    wf = get_wf_gibbs_free_energy(structure, user_kpoints_settings=user_kpoints_settings,
                                  deformations=deformations, vasp_cmd=vasp_cmd, db_file=db_file,
                                  eos=eos, qha_type=qha_type, pressure=pressure, poisson=poisson,
                                  t_min=t_min, t_max=t_max, t_step=t_step, metadata=metadata)

    wf = add_modify_incar(wf, modify_incar_params={"incar_update": {"ENCUT": 600, "EDIFF": 1e-6, "LAECHG": False}})

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


def wf_nudged_elastic_band(structures, parent, c=None):
    """
    Nudged elastic band workflow from the given structures and config dict.

    'is_optimized' default False
    'neb_round' default 1

    Notes:
        Different length of Structure list and "is_optimized" are used to determine the workflow:
        1 structure    # The parent structure & two endpoint indexes provided; need relaxation.
                       # The parent structure & two endpoint indexes provided; no need to relax.
        2 structures   # Two endpoints provided; need to relax two endpoints.
                       # Two relaxed endpoints provided; no need to relax two endpoints.
        >=3 structures # All images including two endpoints are provided.

    Args:
        structures ([Structure]):
            1) The parent structure
            2) Two endpoint structures
            3) All images and endpoints
        parent (Structure): parent structure used to get two endpoints.
        c (dict): workflow config dict, basic format:
            {"fireworks": [],  "common_params": {}, "additional_ep_params": {},
             "additional_neb_params": {}}
    Returns:
        Workflow
    """
    if not(isinstance(structures, list) and len(structures) > 0):
        raise ValueError("structures must be a list of Structure!")

    # config initialization
    c = c or {}
    spec = c.get("common_params", {})
    is_optimized = spec.get("is_optimized", False)
    site_indices = spec.get("site_indices")
    mode = len(structures) if len(structures) < 3 else 3

    if c == {}:  # Default config without yaml provided
        neb_round = 1

    else:  # Check config yaml file
        fw_list = [f['fw'] for f in c.get("fireworks")]
        neb_round = len([f for f in fw_list if "NEBFW" in f])
        if len(structures) == 1:
            mode = 1
            if is_optimized:
                assert len(fw_list) == neb_round + 1
            else:
                assert len(fw_list) == neb_round + 2
        elif len(structures) == 2:
            mode = 2
            if is_optimized:
                assert len(fw_list) == neb_round
            else:
                len(fw_list) == neb_round + 1
        else:  # len(structures) >= 3
            mode = 3

    # Get user_incar_settings, user_kpoints_settings & additional_cust_args
    user_incar_settings = [{}] * (neb_round + 2)
    user_kpoints_settings = [{}] * (neb_round + 2)
    additional_cust_args = [{}] * (neb_round + 2)
    for i in range(1, len(c["fireworks"])+1):
        user_incar_settings[-i] = c["fireworks"][-i].get("user_incar_settings", {})
        user_kpoints_settings[-i] = c["fireworks"][-i].get("user_kpoints_settings", {})
        additional_cust_args[-i] = c["fireworks"][-i].get("additional_cust_args", {})

    kwargs = {"user_incar_settings": user_incar_settings,
              "user_kpoints_settings": user_kpoints_settings,
              "additional_cust_args": additional_cust_args}
    # Assign workflow using mode
    if mode == 1:
        wf = get_wf_neb_from_structure(structure=structures[0], site_indices=site_indices,
                                       additional_spec=spec, **kwargs)
    elif mode == 2:
        wf = get_wf_neb_from_endpoints(parent=parent, endpoints=structures,
                                       additional_spec=spec, **kwargs)
    else:  # mode == 3
        wf = get_wf_neb_from_images(parent=parent, images=structures,
                                    additional_spec=spec, **kwargs)
    return wf
