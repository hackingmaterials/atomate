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


# TODO: @computron: Clarify the config dict -computron
# TODO: @computron: Allow default config dict to be loaded from file -computron

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
    for x in range(4):  # everything but BoltzTrap task
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
    # TODO: @kmathew - this is abuse of the config settings. The config settings should just be
    # high-level things (executive stuff, metadata options) and global (e.g.,
    # user_kpoints_settings.grid_density is not a global parameter). Note that they are also
    # typically capitalized. If the user really wants to tune a workflow (e.g. k-mesh) they should
    # not use the preset workflows which are about trusting the defaults. -computron
    user_kpoints_settings = c.get("user_kpoints_settings", {"grid_density": 7000})
    norm_deformations = c.get("norm_deformations", [-0.01, -0.005, 0.005, 0.01])
    shear_deformations = c.get("shear_deformations", [-0.06, -0.03, 0.03, 0.06])
    optimize_structure = c.get("optimize_structure", True)

    wf = get_wf_elastic_constant(structure, vasp_cmd=vasp_cmd,
                                 norm_deformations=norm_deformations,
                                 shear_deformations=shear_deformations,
                                 db_file=db_file, user_kpoints_settings=user_kpoints_settings,
                                 optimize_structure=optimize_structure)
    wf = add_modify_incar(wf, modify_incar_params={"incar_update":{"ENCUT": 700,
                                                                   "EDIFF": 1e-6,
                                                                   "LAECHG":False}})

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
    # TODO: @kmathew - See my previous comment about config dict params  -computron
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

    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    db_file = c.get("DB_FILE", DB_FILE)

    # TODO: @kmathew - See my previous comment about config dict params  -computron
    user_kpoints_settings = c.get("user_kpoints_settings", {"grid_density": 7000})

    # 21 deformed structures
    deformations = c.get("deformations", [(np.identity(3)*(1+x)).tolist()
                                          for x in np.linspace(-0.1, 0.1, 21)])

    eos = c.get("eos", "vinet")
    qha_type = c.get("qha_type", "debye_model")

    # TODO: @kmathew - See my previous comment about config dict params. Some of these might
    # reasonably be in the config dict, e.g. GIBBS.T_MIN, GIBBS.TMAX  -computron

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

    # TODO: @kmathew - See my previous comment about config dict params.
    # TODO: @kmathew - see also CAPITALIZATION of config_dict params. I fixed some but
    # left some for you. -computron
    c = c or {}
    eos = c.get("eos", "vinet")
    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    db_file = c.get("DB_FILE", DB_FILE)
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
    # TODO: @kmathew - See my previous comment about config dict params.
    # TODO: @kmathew - see also CAPITALIZATION of config_dict params. I fixed some but
    # left some for you. -computron
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
    Nudged elastic band (NEB) workflow from the given structures and config dict.

    'is_optimized' default False
    'neb_round' default 1

    Notes:
        Length of Structure list and "is_optimized" are used to determine the workflow:
        1 structure    # The parent structure & two endpoint indexes provided; need relaxation.
                       # The parent structure & two endpoint indexes provided; no need to relax.
        2 structures   # Two endpoints provided; need to relax two endpoints.
                       # Two relaxed endpoints provided; no need to relax two endpoints.
        >=3 structures # All images including two endpoints are provided.

    Args:
        structures ([Structure]):
            1) The parent structure
            2) Two endpoint structures
            3) An initial NEB path that comprises both images and endpoints
        parent (Structure): parent structure used to get two endpoints.
        c (dict): workflow config dict, basic format:
            {"fireworks": [],  "common_params": {}, "additional_ep_params": {},
             "additional_neb_params": {}}. When the length of structures is 1, "site_indices" key
             must be included in c. Note that "fireworks" is a list corresponding to the order of
             execution.
    Returns:
        Workflow
    """
    if not(isinstance(structures, list) and len(structures) > 0):
        raise ValueError("structures must be a list of Structure objects!")

    # config initialization
    c = c or {}
    # TODO: @shyuep - config dict params are typically capitalized and intended to be global
    # params unless namespaced. e.g. NEB.COMMON_PARAMS if it's only settings for NEB
    # workflows. -computron
    spec = c.get("common_params", {})
    is_optimized = spec.get("is_optimized", False)

    # TODO: @shyuep: the whole point of preset workflows is that they are supposed to be simple.
    # e.g., give a structure or a list of structures, and forget about the rest. Here it looks like
    # one needs to construct some kind of complicated configuration dictionary. Pretty sure no one
    # apart from the people in your group have any idea how to use these functions or set up this
    # config dict (which should not be necessary in the first place). -computron

    if c.get("fireworks"):  # Check config dict file
        fw_list = [f['fw'] for f in c.get("fireworks")]
        neb_round = len([f for f in fw_list if "NEBFW" in f])

        if neb_round < 1:
            raise ValueError("At least one NEB Fireworks (NEBFW) is needed in the config dict!")

        if len(structures) == 1:
            assert "site_indices" in spec, "Site indices not provided in config dict!"
            assert len(spec["site_indices"]) == 2, "Only two site indices should be provided!"

            if is_optimized:
                assert len(fw_list) == neb_round + 1
            else:
                assert len(fw_list) == neb_round + 2
        elif len(structures) == 2:
            if is_optimized:
                assert len(fw_list) == neb_round
            else:
                assert len(fw_list) == neb_round + 1
    else:  # Default settings if config dict is not provided.
        neb_round = 1

    # Get user_incar_settings, user_kpoints_settings & additional_cust_args
    user_incar_settings = [{}] * (neb_round + 2)
    user_kpoints_settings = [{"grid_density": 1000}] * (neb_round + 2)
    additional_cust_args = [{}] * (neb_round + 2)

    # TODO: @shyuep - I have no idea what's going on here -computron
    if "fireworks" in c:
        for i in range(1, len(c["fireworks"])+1):
            user_incar_settings[-i] = c["fireworks"][-i].get("user_incar_settings", {})
            user_kpoints_settings[-i] = c["fireworks"][-i].get("user_kpoints_settings",
                                                               {"grid_density": 1000})
            additional_cust_args[-i] = c["fireworks"][-i].get("additional_cust_args", {})

    kwargs = {"user_incar_settings": user_incar_settings,
              "user_kpoints_settings": user_kpoints_settings,
              "additional_cust_args": additional_cust_args}

    # Assign workflow using the number of given structures
    if len(structures) == 1:
        wf = get_wf_neb_from_structure(structure=structures[0],
                                       additional_spec=spec, **kwargs)
    elif len(structures) == 2:
        wf = get_wf_neb_from_endpoints(parent=parent, endpoints=structures,
                                       additional_spec=spec, **kwargs)
    else:  # len(structures) >= 3
        wf = get_wf_neb_from_images(parent=parent, images=structures,
                                    additional_spec=spec, **kwargs)

    return wf
