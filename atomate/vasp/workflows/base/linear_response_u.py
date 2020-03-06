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
from pymatgen.io.vasp.sets import MPRelaxSet, LinearResponseUSet
from pymatgen.core import Lattice, Structure

import numpy as np

__author__ = ""
__maintainer__ = ""
__email__ = ""
__status__ = "Production"
__date__ = "February 2020"

__linear_response_u_wf_version__ = 0.0

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

def get_wf_linear_response_u(structure, applied_potential_values=[0.0], ground_state_dir=None, c=None, vis=None):
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

    structure = structure.get_primitive_structure(use_site_props=True)

    # using a uuid for book-keeping,
    # in a similar way to other workflows
    uuid = str(uuid4())

    c_defaults = {"vasp_cmd": VASP_CMD, "db_file": DB_FILE}
    if c:
        c.update(c_defaults)
    else:
        c = c_defaults

    # Calculate groundstate
    vis_d = vis.as_dict()
    uis_gs = {"LDAU":False, "LMAXMIX":4, "LORBIT": 11, "ISPIN": 2}
    vis_d.update({"user_incar_settings": uis_gs})
    vis_gs = vis.__class__.from_dict(vis_d)

    if ground_state_dir:
        fws = []
    else:
        fws = [StaticFW(structure=structure, name="initial static", vasp_input_set=vis_gs,
                        vasp_cmd=VASP_CMD, db_file=DB_FILE)]

    uis_ldau = uis_gs.copy()
    uis_ldau.update({"LDAU":True, "LDAUTYPE":3, "LDAUPRINT":2,"ISTART":1})

    if 0.0 in applied_potential_values:
        applied_potential_values = list(applied_potential_values)
        applied_potential_values.pop(applied_potential_values.index(0.0))
        applied_potential_values = np.array(applied_potential_values)

    for v in applied_potential_values:

        sign = 'neg' if str(v)[0] == '-' else 'pos'
        vis_ldau = LinearResponseUSet.from_dict(vis_d)

        for k in ["LDAUL", "LDAUU", "LDAUJ"]:
            val_dict = {}
            if (k == "LDAUL"):
                # for LDAUL
                val_dict.update({"perturb":2})
                for s in vis_ldau.poscar.site_symbols:
                    val_dict.update({s:-1})
                    uis_ldau.update({k:val_dict})
            else:
                # for LDAUU and LDAUJ
                val_dict.update({"perturb":v})
                for s in vis_ldau.poscar.site_symbols:
                    val_dict.update({s:0})
                    uis_ldau.update({k:val_dict})

        # NSCF
        uis_ldau.update({"ICHARG":11})

        vis_d = vis_ldau.as_dict()
        vis_d.update({"user_incar_settings": uis_ldau})
        vis_ldau = LinearResponseUSet.from_dict(vis_d)

        if ground_state_dir:
            parents = []
        else:
            parents=fws[0]
            
        fw = LinearResponseUFW(structure=structure, parents=parents,
                               name="nscf_u_eq_{}{}".format(sign, abs(round(v,6))),
                               vasp_input_set=vis_ldau,
                               vasp_input_set_params = {"user_incar_settings":uis_ldau},
                               additional_files=["WAVECAR","CHGCAR"],
                               prev_calc_dir=ground_state_dir,
                               vasp_cmd=VASP_CMD, db_file=DB_FILE)

        fws.append(fw)

        # SCF
        uis_ldau.update({"ICHARG":0})

        vis_d = vis_ldau.as_dict()
        vis_d.update({"user_incar_settings": uis_ldau})
        vis_ldau = LinearResponseUSet.from_dict(vis_d)

        # NOTE: More efficient to reuse WAVECAR or remove dependency of SCF on NSCF?

        if ground_state_dir:
            parents = []
        else:
            parents=fws[0]
            # parents=fws[-1]
        
        fw = LinearResponseUFW(structure=structure, parents=parents,
                               name="scf_u_eq_{}{}".format(sign, abs(round(v,6))),
                               vasp_input_set=vis_ldau,
                               vasp_input_set_params = {"user_incar_settings":uis_ldau},
                               additional_files=["WAVECAR"],
                               prev_calc_dir=ground_state_dir,
                               vasp_cmd=VASP_CMD, db_file=DB_FILE)
        fws.append(fw)


    fw_analysis = Firework(
        LinearResponseUToDb(
            db_file=DB_FILE, wf_uuid=uuid
        ),
        name="LinearResponseUToDb",
    )

    fws.append(fw_analysis)

    wf = Workflow(fws)
    wf = add_common_powerups(wf, c)

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
