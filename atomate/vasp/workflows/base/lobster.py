#!/usr/bin/env python
# coding: utf-8

"""
This module defines the VASP/Lobster workflows
"""

import itertools
import logging
import os
from collections import OrderedDict
from typing import List

from atomate.common.firetasks.glue_tasks import DeleteFilesPrevFolder
from atomate.vasp.config import VASP_CMD, DB_FILE
from atomate.vasp.fireworks import StaticFW
from atomate.vasp.fireworks.lobster import LobsterFW
from fireworks import Firework
from fireworks.core.firework import Workflow
from pymatgen.core.structure import Structure
from pymatgen.io.lobster import Lobsterin
from pymatgen.io.vasp.sets import LobsterSet

__author__ = "Janine George, Guido Petretto"
__email__ = 'janine.george@uclouvain.be, guido.petretto@uclouvain.be'

LOBSTER_CMD = ">>lobster_cmd<<"
MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
logger = logging.getLogger(__name__)


def get_wf_lobster(structure: Structure, calculationtype: str = 'standard', delete_all_wavecars: bool = True,
                   user_lobsterin_settings: dict = None, user_incar_settings: dict = None,
                   user_kpoints_settings: dict = None, user_supplied_basis: dict = None,
                   isym: int = -1, c: dict = None, additional_outputs: List[str] = None) -> Workflow:
    """
    Creates a workflow for a static vasp calculation followed by a Lobster calculation.

    Args:
        structure: Structure object
        calculationtype: type of the Lobster calculation
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
    staticfw = StaticFW(structure=structure,
                        vasp_input_set=LobsterSet(structure, user_incar_settings=user_incar_settings,
                                                  user_kpoints_settings=user_kpoints_settings, isym=isym),
                        vasp_cmd=vasp_cmd, db_file=db_file)
    fws.append(staticfw)
    fws.append(LobsterFW(structure=structure, parents=staticfw, calculationtype=calculationtype,
                         delete_wavecar=delete_all_wavecars,
                         delete_wavecar_previous_fw=delete_all_wavecars, lobster_cmd=lobster_cmd,
                         db_file=db_file, lobsterin_key_dict=user_lobsterin_settings,
                         user_supplied_basis=user_supplied_basis, handler_group="default",
                         validator_group="strict",
                         additional_outputs=additional_outputs))

    workflow = Workflow(fws, name="LobsterWorkflow")
    return workflow


# TODO: include address_min_basis and address_max_basis correctly
# extend test
def get_wf_lobster_test_basis(structure: Structure, calculationtype: str = 'standard', delete_all_wavecars: bool = True,
                              c: dict = None,
                              address_max_basis: str = os.path.join(MODULE_DIR, "basis_lobster/PBE_54_max_basis.yaml"),
                              address_min_basis: str = os.path.join(MODULE_DIR, "basis_lobster/PBE_54_min_basis.yaml"),
                              user_lobsterin_settings: dict = None,
                              user_incar_settings: dict = None,
                              user_kpoints_settings: dict = None,
                              isym: int = -1, additional_outputs: List[str] = None) -> Workflow:
    """
    creates workflow where all possible basis functions for one compound are tested
    at the end, the user has to decide which projection worked best (e.g., based on chargespilling)
    this is the recommended workflow at the moment!
    Args:
        structure (Structure): structure object that will be used during the run
        calculationtype (str): only "standard" is implemented so far
        delete_all_wavecars (bool): all wavecars wil be delted if True
        c (dict): specifications for wf, e.g. VASP_CMD, LOBSTER_CMD etc.
        address_max_basis (str): address to yaml file including maximum basis set
        address_min_basis (str): address to yaml file including minimum basis set
        user_lobsterin_settings (dict): change lobsterin settings here
        user_incar_settings (dict): change incar settings with this dict
        user_kpoints_settings (dict): change kpoint settings with this dict
        isym (int): isym setting during the VASP calculation, currently lobster can only deal with isym=-1,
            newer versions will be able to deal with 0
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
    inputset = LobsterSet(structure, address_basis_file=address_max_basis,
                          user_incar_settings=user_incar_settings,
                          user_kpoints_settings=user_kpoints_settings, isym=isym)
    # get the basis from dict_max_basis
    potcar_symbols = inputset.potcar_symbols
    max_basis = Lobsterin._get_basis(structure=structure,
                                     potcar_symbols=potcar_symbols, address_basis_file=address_max_basis)
    min_basis = Lobsterin._get_basis(structure=structure,
                                     potcar_symbols=potcar_symbols, address_basis_file=address_min_basis)

    # get information about how many lobster calculations have to be performed to test all possible basis functions
    all_basis = get_all_possible_basis_combinations(min_basis=min_basis, max_basis=max_basis)

    staticfw = StaticFW(structure=structure,
                        vasp_input_set=inputset,
                        vasp_cmd=vasp_cmd, db_file=db_file, name="static")
    fws.append(staticfw)

    # append all lobster calculations that need to be done
    fws_lobster = []
    for ibasis, basis in enumerate(all_basis):
        basis_dict = {}

        for iel, elba in enumerate(basis):
            basplit = elba.split()
            basis_dict[basplit[0]] = " ".join(basplit[1:])

        fws_lobster.append(
            LobsterFW(structure=structure, parents=staticfw, calculationtype=calculationtype,
                      delete_wavecar=delete_all_wavecars,
                      delete_wavecar_previous_fw=False, lobster_cmd=lobster_cmd,
                      db_file=db_file, user_supplied_basis=basis_dict,
                      lobsterin_key_dict=user_lobsterin_settings,
                      handler_group="default", validator_group="strict",
                      name="lobster_calculation_{}".format(ibasis),
                      lobstertodb_kwargs={"additional_fields": {"basis_id": ibasis,
                                                                "number_lobster_runs": len(all_basis)}},
                      additional_outputs=additional_outputs))

    fws.extend(fws_lobster)
    # wavecar from static run is deleted, WAVECARs without symmetry can be huge!
    if delete_all_wavecars:
        t = DeleteFilesPrevFolder(files=["WAVECAR", "WAVECAR.gz"], calc_loc="static")
        final_fw = Firework([t], parents=fws_lobster, name="DelteWavecar")
        fws.append(final_fw)

    workflow = Workflow(fws, name="LobsterWorkflow")
    return workflow


def get_all_possible_basis_combinations(min_basis: list, max_basis: list) -> list:
    """

    Args:
        min_basis: list of basis entries: e.g., ['Si 3p 3s ']
        max_basis: list of basis entries: e.g., ['Si 3p 3s ']

    Returns: all possible combinations of basis functions, e.g. [['Si 3p 3s']]

    """
    max_basis_lists = [x.split() for x in max_basis]
    min_basis_lists = [x.split() for x in min_basis]

    # get all possible basis functions
    basis_dict = OrderedDict({})
    for iel, el in enumerate(max_basis_lists):
        basis_dict[el[0]] = {"fixed": [], "variable": [], "combinations": []}
        for basis in el[1:]:
            if basis in min_basis_lists[iel]:
                basis_dict[el[0]]["fixed"].append(basis)
            if basis not in min_basis_lists[iel]:
                basis_dict[el[0]]["variable"].append(basis)
        for L in range(0, len(basis_dict[el[0]]['variable']) + 1):
            for subset in itertools.combinations(basis_dict[el[0]]['variable'], L):
                basis_dict[el[0]]["combinations"].append(' '.join([el[0]] + basis_dict[el[0]]['fixed'] + list(subset)))

    list_basis = []
    for el, item in basis_dict.items():
        list_basis.append(item['combinations'])

    # get all combinations
    start_basis = list_basis[0]
    if len(list_basis) > 1:
        for iel, el in enumerate(list_basis[1:], 1):
            new_start_basis = []
            for ielbasis, elbasis in enumerate(start_basis):

                for ielbasis2, elbasis2 in enumerate(list_basis[iel]):
                    if type(elbasis) != list:
                        new_start_basis.append([elbasis, elbasis2])
                    else:
                        new_start_basis.append(elbasis.copy() + [elbasis2])
            start_basis = new_start_basis
        return start_basis
    else:
        return [[basis] for basis in start_basis]
