# coding: utf-8

import os

from atomate.vasp.fireworks.core import StaticFW, LinearResponseUFW
from fireworks import Workflow, Firework
from atomate.vasp.powerups import (
    add_tags,
    add_additional_fields_to_taskdocs,
    add_wf_metadata,
    add_common_powerups,
)
from atomate.vasp.workflows.base.core import get_wf

from atomate.vasp.firetasks.parse_outputs import (
    LinearResponseUToDb
)

from atomate.utils.utils import get_logger

logger = get_logger(__name__)

from atomate.vasp.config import VASP_CMD, DB_FILE, ADD_WF_METADATA

from atomate.vasp.workflows.presets.scan import wf_scan_opt
from uuid import uuid4
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet, LinearResponseUSet
from pymatgen.core import Lattice, Structure

import numpy as np

__author__ = ""
__maintainer__ = ""
__email__ = ""
__status__ = "Production"
__date__ = "February 2020"

__linear_response_u_wf_version__ = 0.0

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

def get_wf_linear_response_u(structure, applied_potential_values=[0.0],
                             use_default_uvals=False, ground_state_ldau=True, ground_state_dir=None,
                             c=None, vis=None):
    """
    Args:
        structure: 
        c: Workflow config dict, in the same format
    as in presets/core.py and elsewhere in atomate
        vis: A VaspInputSet to use for the first FW

    Returns: Workflow
    """

    if not structure.is_ordered:
        raise ValueError(
            "Please obtain an ordered approximation of the input structure."
        )

    # using a uuid for book-keeping,
    # in a similar way to other workflows
    uuid = str(uuid4())

    c_defaults = {"vasp_cmd": VASP_CMD, "db_file": DB_FILE}
    if c:
        c.update(c_defaults)
    else:
        c = c_defaults

    # Calculate groundstate

    # ground state user incar settings
    uis_gs = {"LDAU":False, "LMAXMIX":4, "LORBIT": 11, "ISPIN": 2}

    uis_ldau = uis_gs.copy()
    uis_ldau.update({"LDAU":True, "LDAUTYPE":3, "LDAUPRINT":2})

    # Load default U values
    vis_params = {"user_incar_settings":{}}
    set_default_uvals = MPRelaxSet(structure=structure, sort_structure=False, **vis_params)
    incar_dict_default_u = set_default_uvals.incar.as_dict()
    sitesym_default_u = set_default_uvals.poscar.site_symbols
    default_uvals = {}
    if 'LDAUU' in incar_dict_default_u.keys():
        uvals = incar_dict_default_u['LDAUU']
        lvals = incar_dict_default_u['LDAUL']
        for sym, u, l in zip(sitesym_default_u, uvals, lvals):
            default_uvals.update({sym:{'LDAUU':u, 'LDAUL':l}})

    # initialize vasp input set
    vis_params = {"user_incar_settings": uis_ldau.copy()}
    vis_ldau = LinearResponseUSet(structure=structure, **vis_params)

    for k in ["LDAUL", "LDAUU", "LDAUJ"]:
        val_dict = {}
        if (k == "LDAUL"):
            # for LDAUL
            val_dict.update({"perturb":2})             # FIXME: shouldn't hard code LDAUL
            for s in vis_ldau.poscar.site_symbols:
                l = -1
                if use_default_uvals:
                    if s in default_uvals.keys():
                        if k in default_uvals[s].keys():
                            l = default_uvals[s][k]
                val_dict.update({s:l})
            uis_ldau.update({k:val_dict})
        else:
            # for LDAUU and LDAUJ
            val_dict.update({"perturb":0})
            for s in vis_ldau.poscar.site_symbols:
                v = 0
                if use_default_uvals:
                    if s in default_uvals.keys():
                        if 'LDAUU' in default_uvals[s].keys():
                            v = default_uvals[s]['LDAUU']
                val_dict.update({s:v})
            uis_ldau.update({k:val_dict})

    if ground_state_ldau:
        uis_gs = uis_ldau.copy()
        vis_params = {"user_incar_settings": uis_gs.copy()}
        vis_gs = LinearResponseUSet(structure=structure, **vis_params)
        fw_gs = LinearResponseUFW(structure=structure, name="initial static", vasp_input_set=vis_gs,
                                  vasp_cmd=VASP_CMD, db_file=DB_FILE)
    else:
        vis_params = {"user_incar_settings": uis_gs.copy()}
        vis_gs = MPStaticSet(structure=structure, sort_structure=False, **vis_params)
        fw_gs = StaticFW(structure=structure, name="initial static", vasp_input_set=vis_gs,
                         vasp_cmd=VASP_CMD, db_file=DB_FILE)

    if ground_state_dir:
        fws = []
    else:
        fws = [fw_gs]

    if 0.0 in applied_potential_values:
        applied_potential_values = list(applied_potential_values)
        applied_potential_values.pop(applied_potential_values.index(0.0))
        applied_potential_values = np.array(applied_potential_values)

    for v in applied_potential_values:

        sign = 'neg' if str(v)[0] == '-' else 'pos'

        # Update applied potential to U and J
        uis_ldau.update({"ISTART":1})

        # Update perturbation potential for U and J
        for k in ["LDAUU", "LDAUJ"]:
            # for LDAUU and LDAUJ
            val_dict.update({"perturb":v})
            uis_ldau.update({k:val_dict})

        vis_params = {"user_incar_settings": uis_ldau.copy()}
        vis_ldau = LinearResponseUSet(structure=structure, **vis_params)

        # NSCF
        uis_ldau.update({"ICHARG":11})

        vis_params = {"user_incar_settings": uis_ldau.copy()}
        vis_ldau = LinearResponseUSet(structure=structure, **vis_params)

        if ground_state_dir:
            parents = []
        else:
            parents=fws[0]

        fw = LinearResponseUFW(structure=structure, parents=parents,
                               name="nscf_u_eq_{}{}".format(sign, abs(round(v,6))),
                               vasp_input_set=vis_ldau,
                               additional_files=["WAVECAR","CHGCAR"],
                               prev_calc_dir=ground_state_dir,
                               vasp_cmd=VASP_CMD, db_file=DB_FILE)

        fws.append(fw)

        # SCF
        uis_ldau.update({"ICHARG":0})

        vis_params = {"user_incar_settings": uis_ldau.copy()}
        vis_ldau = LinearResponseUSet(structure=structure, **vis_params)

        # NOTE: More efficient to reuse WAVECAR or remove dependency of SCF on NSCF?

        if ground_state_dir:
            parents = []
        else:
            parents=fws[0]
            # parents=fws[-1]

        fw = LinearResponseUFW(structure=structure, parents=parents,
                               name="scf_u_eq_{}{}".format(sign, abs(round(v,6))),
                               vasp_input_set=vis_ldau,
                               additional_files=["WAVECAR"],
                               prev_calc_dir=ground_state_dir,
                               vasp_cmd=VASP_CMD, db_file=DB_FILE)
        fws.append(fw)

    wf = Workflow(fws)
    # Needed here?
    wf = add_common_powerups(wf, c)

    # fw_analysis = Firework(
    #     LinearResponseUToDb(
    #         db_file=DB_FILE, wf_uuid=uuid
    #     ),
    #     name="LinearResponseUToDb",
    # )

    # wf.append_wf(Workflow.from_Firework(fw_analysis), wf.leaf_fw_ids)

    # wf = add_common_powerups(wf, c)

    # if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
    #     wf = add_wf_metadata(wf, structure)

    wf = add_additional_fields_to_taskdocs(
        wf,
        {
            "wf_meta": {
                "wf_uuid": uuid,
                "wf_name": "linear_response_u",
                "wf_version": __linear_response_u_wf_version__,
            }
        },
    )

    return wf
