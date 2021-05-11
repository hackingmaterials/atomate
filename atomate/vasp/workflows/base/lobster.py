#!/usr/bin/env python

"""
This module defines the VASP/Lobster workflows
"""

import logging
import os
from typing import List, Optional

from fireworks import Firework
from fireworks.core.firework import Workflow

from atomate.common.firetasks.glue_tasks import DeleteFilesPrevFolder
from atomate.vasp.config import VASP_CMD, DB_FILE, LOBSTER_CMD
from atomate.vasp.fireworks import StaticFW
from atomate.vasp.fireworks.lobster import LobsterFW
from pymatgen.core.structure import Structure
from pymatgen.io.lobster import Lobsterin
from pymatgen.io.vasp.sets import LobsterSet

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
) -> Workflow:
    """
    Creates a workflow for a static Vasp calculation followed by a Lobster calculation.

    Args:
        structure: Structure object
        calculation_type: type of the Lobster calculation
        delete_all_wavecars: if True, all WAVECARs are deleted
        user_lobsterin_settings (dict): dict to set additional lobsterin settings
        user_incar_settings (dict): dict to set additional things in INCAR
        user_kpoints_settings (dict): dict to set additional things in KPOINTS
        user_supplied_basis (dict): dict to supply basis functions for each element type
        isym (int): isym setting during the vasp calculation, currently lobster can only deal with isym=-1
        c (dict): configurations dict which can include "VASP_CMD", "LOBSTER_CMD", "DB_FILE"
        additional_outputs (list): list of additional files to be stored in the
            results DB. They will be stored as files in gridfs. Examples are:
            "ICOHPLIST.lobster" or "DOSCAR.lobster". Note that the file name
            should be given with the full name and the correct capitalization.

    Returns: Workflow


    """

    c = c or {}

    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    lobster_cmd = c.get("LOBSTER_CMD", LOBSTER_CMD)
    db_file = c.get("DB_FILE", DB_FILE)

    fws = []
    staticfw = StaticFW(
        structure=structure,
        vasp_input_set=LobsterSet(
            structure,
            user_incar_settings=user_incar_settings,
            user_kpoints_settings=user_kpoints_settings,
            isym=isym,
        ),
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
) -> Workflow:
    """
    creates workflow where all possible basis functions for one compound are tested
    at the end, the user has to decide which projection worked best (e.g., based on chargespilling)
    this is the recommended workflow at the moment!
    Args:
        structure (Structure): structure object that will be used during the run
        calculation_type (str): only "standard" is implemented so far
        delete_all_wavecars (bool): all wavecars wil be deleted if True
        c (dict): specifications for wf, e.g. VASP_CMD, LOBSTER_CMD etc.
        address_max_basis (str): address to yaml file including maximum basis set (otherwise predefined file)
        address_min_basis (str): address to yaml file including minimum basis set (otherwise predefined file)
        user_lobsterin_settings (dict): change lobsterin settings here
        user_incar_settings (dict): change incar settings with this dict
        user_kpoints_settings (dict): change kpoint settings with this dict
        isym (int): isym setting during the VASP calculation, currently lobster can only deal with isym=-1 and isym=0
        additional_outputs (list): list of additional files to be stored in the
            results DB. They will be stored as files in gridfs. Examples are:
            "ICOHPLIST.lobster" or "DOSCAR.lobster". Note that the file name
            should be given with the full name and the correct capitalization.
    Returns:
    """

    c = c or {}

    vasp_cmd = c.get("VASP_CMD", VASP_CMD)
    lobster_cmd = c.get("LOBSTER_CMD", LOBSTER_CMD)
    db_file = c.get("DB_FILE", DB_FILE)

    fws = []
    # get the relevant potcar files!
    inputset = LobsterSet(
        structure,
        address_basis_file=address_max_basis,
        user_incar_settings=user_incar_settings,
        user_kpoints_settings=user_kpoints_settings,
        isym=isym,
    )
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
                validator_group="strict",
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
