#!/usr/bin/env python

"""
This module defines the VASP/Lobster workflows
"""

import logging
import os
from typing import List, Optional

from fireworks import Firework, Workflow
from pymatgen.core import Structure
from pymatgen.io.lobster import Lobsterin
from pymatgen.io.vasp.sets import LobsterSet

from atomate.common.firetasks.glue_tasks import DeleteFilesPrevFolder
from atomate.vasp.config import DB_FILE, LOBSTER_CMD, VASP_CMD
from atomate.vasp.fireworks import OptimizeFW, StaticFW
from atomate.vasp.fireworks.lobster import LobsterFW

__author__ = "Janine George, Guido Petretto"
__email__ = "janine.george@uclouvain.be, guido.petretto@uclouvain.be"

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
logger = logging.getLogger(__name__)


def get_wf_lobster(
    structure: Structure,
    calculation_type: str = "standard",
    delete_all_wavecars: bool = True,
    user_lobsterin_settings: dict = None,
    user_incar_settings: dict = None,
    user_kpoints_settings: dict = None,
    user_supplied_basis: dict = None,
    isym: int = 0,
    c: dict = None,
    additional_outputs: List[str] = None,
    additional_optimization: bool = False,
    user_incar_settings_optimization: dict = None,
    user_kpoints_settings_optimization: dict = None,
) -> Workflow:
    """
    Creates a workflow for a static Vasp calculation followed by a Lobster calculation.

    Args:
        structure: Structure object
        calculation_type: type of the Lobster calculation
        delete_all_wavecars: if True, all WAVECARs are deleted
        user_lobsterin_settings (dict): dict to set additional lobsterin settings
        user_incar_settings (dict): dict to set additional things in INCAR, only for LobsterCalc not optimization
        user_kpoints_settings (dict): dict to set additional things in KPOINTS, only for LobsterCalc not optimization
        user_supplied_basis (dict): dict to supply basis functions for each element type
        isym (int): isym setting during the vasp calculation, currently lobster can only deal with isym=-1
        c (dict): configurations dict which can include "VASP_CMD", "LOBSTER_CMD", "DB_FILE"
        additional_outputs (list): list of additional files to be stored in the
            results DB. They will be stored as files in gridfs. Examples are:
            "ICOHPLIST.lobster" or "DOSCAR.lobster". Note that the file name
            should be given with the full name and the correct capitalization.
        additional_optimization (bool): determines if an optimization is performed
        user_incar_settings_optimization (dict): change incar settin with this dict for optimization
        user_kpoints_settings_optimization (dict): change incar settin with this dict for optimization
    Returns: Workflow


    """

    c = c or {}

    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    lobster_cmd = c.get("LOBSTER_CMD", LOBSTER_CMD)
    db_file = c.get("DB_FILE", DB_FILE)

    fws = []

    lobster_set = LobsterSet(
        structure,
        user_incar_settings=user_incar_settings,
        user_kpoints_settings=user_kpoints_settings,
        isym=isym,
    )
    if additional_optimization:

        # add an additional optimization firework
        optmize_fw = OptimizeFW(
            structure,
            override_default_vasp_params={
                "user_potcar_functional": "PBE_54",
                "user_potcar_settings": {"W": "W_sv"},
                "user_kpoints_settings": user_kpoints_settings_optimization,
                "user_incar_settings": user_incar_settings_optimization,
            },
            vasp_cmd=vasp_cmd,
            db_file=db_file,
        )
        fws.append(optmize_fw)
        user_incar_settings_here = lobster_set.incar.as_dict()
        user_incar_settings_here.update(
            {"ISTART": None, "LAECHG": None, "LCHARG": None, "LVHAR": None}
        )
        if user_incar_settings is not None:
            user_incar_settings_here.update(user_incar_settings)
        user_kpoints_settings_here = lobster_set.kpoints.as_dict()
        if user_kpoints_settings is not None:
            user_incar_settings_here.update(user_kpoints_settings)
        staticfw = StaticFW(
            structure=structure,
            vasp_input_set=lobster_set,
            vasp_input_set_params={
                "user_incar_settings": user_incar_settings_here,
                "user_potcar_functional": "PBE_54",
                "user_potcar_settings": {"W": "W_sv"},
                "user_kpoints_settings": user_kpoints_settings_here,
            },
            vasp_cmd=vasp_cmd,
            db_file=db_file,
            parents=optmize_fw,
        )

    else:

        staticfw = StaticFW(
            structure=structure,
            vasp_input_set=lobster_set,
            vasp_cmd=vasp_cmd,
            db_file=db_file,
        )
    fws.append(staticfw)
    fws.append(
        LobsterFW(
            structure=structure,
            parents=staticfw,
            calculation_type=calculation_type,
            delete_wavecar=delete_all_wavecars,
            delete_wavecar_previous_fw=delete_all_wavecars,
            lobster_cmd=lobster_cmd,
            db_file=db_file,
            lobsterin_key_dict=user_lobsterin_settings,
            user_supplied_basis=user_supplied_basis,
            handler_group="default",
            validator_group="strict",
            additional_outputs=additional_outputs,
        )
    )

    workflow = Workflow(fws, name="LobsterWorkflow")
    return workflow


def get_wf_lobster_test_basis(
    structure: Structure,
    calculation_type: str = "standard",
    delete_all_wavecars: bool = True,
    c: dict = None,
    address_max_basis: Optional[str] = None,
    address_min_basis: Optional[str] = None,
    user_lobsterin_settings: dict = None,
    user_incar_settings: dict = None,
    user_kpoints_settings: dict = None,
    isym: int = 0,
    additional_outputs: List[str] = None,
    additional_optimization: bool = False,
    user_incar_settings_optimization: dict = None,
    user_kpoints_settings_optimization: dict = None,
) -> Workflow:
    """
    creates workflow where all possible basis functions for one compound are tested
    at the end, the user has to decide which projection worked best (e.g., based on chargespilling)
    this is the recommended workflow at the moment!
    Args:
        structure (Structure): structure object that will be used during the run
        calculation_type (str): only "standard" is implemented so far
        delete_all_wavecars (bool): all wavecars will be deleted if True
        c (dict): specifications for wf, e.g. VASP_CMD, LOBSTER_CMD etc.
        address_max_basis (str): address to yaml file including maximum basis set (otherwise predefined file)
        address_min_basis (str): address to yaml file including minimum basis set (otherwise predefined file)
        user_lobsterin_settings (dict): change lobsterin settings here
        user_incar_settings (dict): change incar settings with this dict, only for Lobster calc and not optimization
        user_kpoints_settings (dict): change KPOINT settings with this dict, only for Lobster calc and not optimization
        isym (int): isym setting during the VASP calculation, currently lobster can only deal with isym=-1 and isym=0
        additional_outputs (list): list of additional files to be stored in the
            results DB. They will be stored as files in gridfs. Examples are:
            "ICOHPLIST.lobster" or "DOSCAR.lobster". Note that the file name
            should be given with the full name and the correct capitalization.
        additional_optimization (bool): determines if an optimization is performed
        user_incar_settings_optimization (dict): change incar settin with this dict for optimization
        user_kpoints_settings_optimization (dict): change incar settin with this dict for optimization
    Returns:
    """

    c = c or {}

    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    lobster_cmd = c.get("LOBSTER_CMD", LOBSTER_CMD)
    db_file = c.get("DB_FILE", DB_FILE)

    # get the relevant potcar files!
    inputset = LobsterSet(
        structure,
        address_basis_file=address_max_basis,
        user_incar_settings=user_incar_settings,
        user_kpoints_settings=user_kpoints_settings,
        isym=isym,
    )

    fws = []
    if additional_optimization:
        optmize_fw = OptimizeFW(
            structure,
            override_default_vasp_params={
                "user_potcar_functional": "PBE_54",
                "user_potcar_settings": {"W": "W_sv"},
                "user_kpoints_settings": user_kpoints_settings_optimization,
                "user_incar_settings": user_incar_settings_optimization,
            },
            vasp_cmd=vasp_cmd,
            db_file=db_file,
        )
        fws.append(optmize_fw)

    # get the basis from dict_max_basis
    potcar_symbols = inputset.potcar_symbols

    # will get all possible basis functions that have to be tested
    if address_max_basis is None and address_min_basis is None:
        list_basis_dict = Lobsterin.get_all_possible_basis_functions(
            structure=structure, potcar_symbols=potcar_symbols
        )
    elif address_max_basis is not None and address_min_basis is None:
        list_basis_dict = Lobsterin.get_all_possible_basis_functions(
            structure=structure,
            potcar_symbols=potcar_symbols,
            address_basis_file_max=address_max_basis,
        )
    elif address_min_basis is not None and address_max_basis is None:
        list_basis_dict = Lobsterin.get_all_possible_basis_functions(
            structure=structure,
            potcar_symbols=potcar_symbols,
            address_basis_file_min=address_min_basis,
        )
    elif address_min_basis is not None and address_max_basis is not None:
        list_basis_dict = Lobsterin.get_all_possible_basis_functions(
            structure=structure,
            potcar_symbols=potcar_symbols,
            address_basis_file_max=address_max_basis,
            address_basis_file_min=address_min_basis,
        )
    if additional_optimization:

        user_incar_settings_here = inputset.incar.as_dict()
        user_incar_settings_here.update(
            {"ISTART": None, "LAECHG": None, "LCHARG": None, "LVHAR": None}
        )
        if user_incar_settings is not None:
            user_incar_settings_here.update(user_incar_settings)
        user_kpoints_settings_here = inputset.kpoints.as_dict()
        if user_kpoints_settings is not None:
            user_kpoints_settings_here.update(user_kpoints_settings)
        staticfw = StaticFW(
            structure=structure,
            vasp_input_set_params={
                "user_incar_settings": user_incar_settings_here,
                "user_potcar_functional": "PBE_54",
                "user_potcar_settings": {"W": "W_sv"},
                "user_kpoints_settings": user_kpoints_settings_here,
            },
            vasp_cmd=vasp_cmd,
            db_file=db_file,
            name="static",
            parents=optmize_fw,
        )
    else:
        staticfw = StaticFW(
            structure=structure,
            vasp_input_set=inputset,
            vasp_cmd=vasp_cmd,
            db_file=db_file,
            name="static",
        )

    fws.append(staticfw)

    # append all lobster calculations that need to be done
    fws_lobster = []
    for ibasis, basis_dict in enumerate(list_basis_dict):
        fws_lobster.append(
            LobsterFW(
                structure=structure,
                parents=staticfw,
                calculation_type=calculation_type,
                delete_wavecar=delete_all_wavecars,
                delete_wavecar_previous_fw=False,
                lobster_cmd=lobster_cmd,
                db_file=db_file,
                user_supplied_basis=basis_dict,
                lobsterin_key_dict=user_lobsterin_settings,
                handler_group="default",
                validator_group="default",
                name=f"lobster_calculation_{ibasis}",
                lobstertodb_kwargs={
                    "basis_id": ibasis,
                    "number_lobster_runs": len(list_basis_dict),
                },
                additional_outputs=additional_outputs,
            )
        )

    fws.extend(fws_lobster)
    # wavecar from static run is deleted, WAVECARs without symmetry can be huge!
    if delete_all_wavecars:
        t = DeleteFilesPrevFolder(files=["WAVECAR", "WAVECAR.gz"], calc_loc="static")
        final_fw = Firework([t], parents=fws_lobster, name="DelteWavecar")
        fws.append(final_fw)

    workflow = Workflow(fws, name="LobsterWorkflow")
    return workflow
