#!/usr/bin/env python
# coding: utf-8

"""
This module defines the Nudged Elastic Band (NEB) workflow.
"""

from datetime import datetime

from pymatgen_diffusion.neb.io import get_endpoints_from_index

from fireworks.core.firework import Workflow

from atomate.vasp.fireworks.core import NEBFW, NEBRelaxationFW

__author__ = "Hanmei Tang, Iek-Heng Chu"
__email__ = 'hat003@eng.ucsd.edu, ihchu@eng.ucsd.edu'


def _update_spec(additional_spec):
    """
    Update spec to overwrite default settings.
    Args:
        additional_spec (dict): user spec settings.
            "is_optimized" (bool): True denotes the provided structures are optimized.
                Otherwise run relaxation for given structures. Default False.
            "interpolation_type" (str): Method to do image interpolation from two endpoints.
                Default "IDPP", otherwise "linear".
            "idpp_species" (str): Species used in IDPP interpolation method.
            "sort_tol" (float): Distance tolerance (in Angstrom) used to match the atomic indices
                between start and end structures. If it is set 0, then no sorting will be performed.
            "d_img" (float): Distance between images, in Angstrom. If "IMAGES" is not provided in
                user_incar_settings, this will be used to calculate the number of images.
                Default 0.7 Angstrom.
            "wf_name" (str): An appropriate and unique name for the workflow.
                The workflow result will be transferred to ">>run_dest_root<</wf_name".
    Returns:
        spec dict
    """
    additional_spec = additional_spec or {}
    default_spec = {"is_optimized": False,
                    "interpolation_type": "IDPP",
                    "idpp_species": None,
                    "sort_tol": 0,
                    "d_img": 0.7,
                    "wf_name": datetime.utcnow().strftime('%Y-%m-%d-%H-%M-%S-%f')}
    default_spec.update(additional_spec)
    return default_spec


def get_wf_neb_from_structure(structure, site_indices, user_incar_settings=None,
                              additional_spec=None, user_kpoints_settings=None,
                              additional_cust_args=None):
    """
    Get the CI-NEB workflow using a parent structure.
    Workflow: (parent relax) -- Endpoints relax -- NEB_1 -- NEB_2 - ... - NEB_r
              parent not optimized: parent--ep--neb(r)
              parent is optimized: ep--neb(r)
    Args:
        structure (Structure): The parent structure.
        site_indices (list): The two vacancy site indices.
        user_incar_settings([dict]): Additional user_incar_settings.
            Always arranged as: "parent", "ep_relax", "neb1", "neb2" etc., which contains at least
            three items. The first dict is for parent structure relaxation, the second dict is for
            endpoints relaxation, and the rest are for NEB calculations. For example, [{}, {},
            {"IOPT": 7}, {"IOPT": 1}]. Besides, user_incar_settings is used to determine how many
            NEB rounds will be, default [{}, {}, {}].
        additional_spec (dict): user spec settings to overwrite default_spec.
        user_kpoints_settings ([dict]): Additional user_kpoints_settings, which contains at at least
            three items. The structure is the same with user_incar_settings. For example, [{}, {},
            {"grid_density": 100}] for the workflow from the parent structure relaxation, then the
            endpoint relaxation followed by one-round NEB simulation. Default value depends on
            chosen VaspInputSet.
        additional_cust_args ([dict]): Optional parameters for RunVaspCustodian, same structure with
            user_incar_settings and user_kpoints_settings.
    Returns:
        Workflow
    """
    spec = _update_spec(additional_spec)
    site_indices = spec["site_indices"]
    is_optimized = spec["is_optimized"]
    wf_name = spec["wf_name"]

    ep0, ep1 = get_endpoints_from_index(structure, site_indices)
    spec["parent"] = structure.as_dict()
    spec["eps"] = [ep0.as_dict(), ep1.as_dict()]

    # Assume one round NEB if user_incar_settings not provided.
    user_incar_settings = user_incar_settings or [{}, {}, {}]
    neb_round = len(user_incar_settings[2:])
    user_kpoints_settings = user_kpoints_settings or [{}] * (neb_round + 2)
    additional_cust_args = additional_cust_args or [{}] * (neb_round + 2)
    for incar in user_incar_settings[2:]:
        if incar.get("IMAGES"):
            # If "incar_images" shows up, the number of images is pre-defined
            spec["incar_images"] = incar["IMAGES"]
            break

    if is_optimized:  # Start from endpoints
        endpoints = get_endpoints_from_index(structure, site_indices)
        endpoints_dict = [e.as_dict() for e in endpoints]
        spec["eps"] = endpoints_dict
        neb_fws, rlx_fws = [], []

        # Get neb fireworks.
        for n in range(neb_round):
            fw = NEBFW(spec=spec, neb_label=str(n + 1), from_images=False,
                       user_incar_settings=user_incar_settings[n + 2],
                       user_kpoints_settings=user_kpoints_settings[n + 2],
                       additional_cust_args=additional_cust_args[n + 2])
            neb_fws.append(fw)
        # Get relax fireworks
        for label in ["ep0", "ep1"]:
            fw = NEBRelaxationFW(spec=spec, label=label,
                                 user_incar_settings=user_incar_settings[1],
                                 user_kpoints_settings=user_kpoints_settings[1],
                                 additional_cust_args=additional_cust_args[1])
            rlx_fws.append(fw)
        # Build fireworks link
        links = {rlx_fws[0]: [neb_fws[0]], rlx_fws[1]: [neb_fws[0]]}

    else:  # Start from perfect structure
        spec["parent"] = structure.as_dict()
        neb_fws, rlx_fws = [], []

        # Get neb fireworks.
        for n in range(neb_round):
            fw = NEBFW(spec=spec, neb_label=str(n + 1), from_images=False,
                       user_incar_settings=user_incar_settings[n + 2],
                       user_kpoints_settings=user_kpoints_settings[n + 2],
                       additional_cust_args=additional_cust_args[n + 2])
            neb_fws.append(fw)
        # Get relaxation fireworks.
        rlx_fws.append(NEBRelaxationFW(spec=spec, label="parent",
                                       user_incar_settings=user_incar_settings[0],
                                       user_kpoints_settings=user_kpoints_settings[0],
                                       additional_cust_args=additional_cust_args[0]))

        for i, label in enumerate(["ep0", "ep1"]):
            fw = NEBRelaxationFW(spec=spec, label=label,
                                 user_incar_settings=user_incar_settings[1],
                                 user_kpoints_settings=user_kpoints_settings[1],
                                 additional_cust_args=additional_cust_args[1])
            rlx_fws.append(fw)

        # Build fireworks link
        links = {rlx_fws[0]: [rlx_fws[1], rlx_fws[2]],
                 rlx_fws[1]: [neb_fws[0]],
                 rlx_fws[2]: [neb_fws[0]]}

    # Put all fireworks together with link
    fws = rlx_fws + neb_fws
    if neb_round >= 2:
        for r in range(1, neb_round):
            links[neb_fws[r - 1]] = [neb_fws[r]]
    workflow = Workflow(fws, links_dict=links, name=wf_name)
    return workflow


def get_wf_neb_from_endpoints(parent, endpoints, user_incar_settings=None, additional_spec=None,
                              user_kpoints_settings=None, additional_cust_args=None):
    """
    Get a CI-NEB workflow from given endpoints.
    Workflow: (Endpoints relax -- ) NEB_1 -- NEB_2 - ... - NEB_r
              endpoints not optimized: ep--neb(r)
              endpoints are optimized: neb(r)

    Args:
        parent (Structure): parent structure.
        endpoints (list[Structure]): The endpoint structures.
        user_incar_settings([dict]): Additional user_incar_settings.
            Always arranged as: "parent", "ep_relax", "neb1", "neb2" etc., which contains at least
            three items. The first dict is for parent structure relaxation, the second dict is for
            endpoints relaxation, and the rest are for NEB calculations. For example, [{}, {},
            {"IOPT": 7}, {"IOPT": 1}]. Besides, user_incar_settings is used to determine how many
            NEB rounds will be, default [{}, {}, {}].
        additional_spec (dict): user spec settings to overwrite default_spec.
        user_kpoints_settings ([dict]): Additional user_kpoints_settings, which contains at at least
            three items. The structure is the same with user_incar_settings. For example, [{}, {},
            {"grid_density": 100}] for the workflow from the parent structure relaxation, then the
            endpoint relaxation followed by one-round NEB simulation. Default value depends on
            chosen VaspInputSet.
        additional_cust_args ([dict]): Optional parameters for RunVaspCustodian, same structure with
            user_incar_settings and user_kpoints_settings.
    Returns:
        Workflow
    """
    spec = _update_spec(additional_spec)
    spec["parent"] = parent.as_dict()
    spec["eps"] = [s.as_dict() for s in endpoints]

    wf_name = spec["wf_name"]
    is_optimized = spec["is_optimized"]

    # Assume one round NEB if user_incar_settings not provided.
    user_incar_settings = user_incar_settings or [{}, {}, {}]
    neb_round = len(user_incar_settings[2:])
    user_kpoints_settings = user_kpoints_settings or [{}] * (neb_round + 2)
    additional_cust_args = additional_cust_args or [{}] * (neb_round + 2)
    for incar in user_incar_settings[2:]:
        if incar.get("IMAGES"):
            # If "incar_images" shows up, the number of images is pre-defined
            spec["incar_images"] = incar["IMAGES"]
            break

    neb_fws = []
    for n in range(neb_round):
        fw = NEBFW(spec=spec, neb_label=str(n + 1), from_images=False,
                   user_incar_settings=user_incar_settings[n + 2],
                   user_kpoints_settings=user_kpoints_settings[n + 2],
                   additional_cust_args=additional_cust_args[n + 2])
        neb_fws.append(fw)

    workflow = Workflow(neb_fws, name=wf_name)

    # Add endpoints relaxation if structures not optimized.
    if not is_optimized:
        ep_fws = []
        for i in ["ep0", "ep1"]:
            fw = NEBRelaxationFW(spec=spec, label=i, user_incar_settings=user_incar_settings[1],
                                 user_kpoints_settings=user_kpoints_settings[1],
                                 additional_cust_args=additional_cust_args[1])
            ep_fws.append(fw)

        # Build fireworks link
        fws = ep_fws + neb_fws
        links = {ep_fws[0]: [neb_fws[0]], ep_fws[1]: [neb_fws[0]]}
        for r in range(1, neb_round):
            links[neb_fws[r - 1]] = [neb_fws[r]]

        workflow = Workflow(fws, links_dict=links, name=wf_name)

    return workflow


def get_wf_neb_from_images(parent, images, user_incar_settings, additional_spec=None,
                           user_kpoints_settings=None, additional_cust_args=None):
    """
    Get a CI-NEB workflow from given images.
    Workflow: NEB_1 -- NEB_2 - ... - NEB_n

    Args:
        parent (Structure): parent structure.
        images ([Structure]): All images and two endpoints.
        user_incar_settings([dict]): Additional user_incar_settings.
            Always arranged as: "parent", "ep_relax", "neb1", "neb2" etc., which contains at least
            three items. The first dict is for parent structure relaxation, the second dict is for
            endpoints relaxation, and the rest are for NEB calculations. For example, [{}, {},
            {"IOPT": 7}, {"IOPT": 1}]. Besides, user_incar_settings is used to determine how many
            NEB rounds will be, default [{}, {}, {}].
        additional_spec (dict): user spec settings to overwrite default_spec.
        user_kpoints_settings ([dict]): Additional user_kpoints_settings, which contains at at least
            three items. The structure is the same with user_incar_settings. For example, [{}, {},
            {"grid_density": 100}] for the workflow from the parent structure relaxation, then the
            endpoint relaxation followed by one-round NEB simulation. Default value depends on
            chosen VaspInputSet.
        additional_cust_args ([dict]): Optional parameters for RunVaspCustodian, same structure with
            user_incar_settings and user_kpoints_settings.

    Returns:
        Workflow
    """
    spec = _update_spec(additional_spec)
    spec["parent"] = parent.as_dict()
    assert isinstance(images, list) and len(images) >= 3
    spec["neb"] = [[s.as_dict() for s in images]]
    spec["_queueadapter"] = {"nnodes": str(len(images) - 2), "nodes": str(len(images) - 2)}
    wf_name = spec["wf_name"]

    # Assume one round NEB if user_incar_settings not provided.
    user_incar_settings = user_incar_settings or [{}, {}, {}]
    neb_round = len(user_incar_settings[2:])
    user_kpoints_settings = user_kpoints_settings or [{}] * (neb_round + 2)
    additional_cust_args = additional_cust_args or [{}] * (neb_round + 2)

    fws = []
    # Get neb fireworks.
    for n in range(neb_round):
        fw = NEBFW(spec=spec, neb_label=str(n + 1), from_images=True,
                   user_incar_settings=user_incar_settings[n + 2],
                   user_kpoints_settings=user_kpoints_settings[n + 2],
                   additional_cust_args=additional_cust_args[n + 2])
        fws.append(fw)

    # Build fireworks link
    links_dict = {}
    if neb_round >= 2:
        for i in range(neb_round - 1):
            links_dict[fws[i]] = [fws[i + 1]]
    workflow = Workflow(fws, name=wf_name, links_dict=links_dict)

    return workflow


if __name__ == "__main__":
    pass
