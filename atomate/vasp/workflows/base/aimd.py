import math
import warnings
from copy import deepcopy
from typing import Dict, List, Optional, Union
import numpy as np
from monty.serialization import loadfn, dumpfn

from atomate.utils.utils import get_logger
from atomate.vasp.config import DB_FILE, VASP_CMD
from atomate.vasp.firetasks import pass_vasp_result
from atomate.vasp.analysis.lattice_dynamics import (
    FIT_METHOD,
    MESH_DENSITY,
)
from atomate.vasp.fireworks.core import TransmuterFW, OptimizeFW
from atomate.vasp.fireworks.aimd import ARCMDFW, CollectMDSegmentsFW
from atomate.common.powerups import add_additional_fields_to_taskdocs
from fireworks import Workflow
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.sets import ARCMDSet, VaspInputSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.advanced_transformations import (
    CubicSupercellTransformation,
)

__author__ = "Junsoo Park"
__email__ = "junsoo.park@nasa.gov"
__date__ = "December 2023"

logger = get_logger(__name__)

vasp_to_db_params = {
    "store_volumetric_data": tuple(),
    "vasp_drone_params": {"parse_bader": False, "parse_locpot": False}
}

_DEFAULT_SETTINGS = {"VASP_CMD": VASP_CMD, "DB_FILE": DB_FILE}
_WF_VERSION = 0.1


def get_aimd_wf(
    structure: Structure,
    common_settings: Dict = None,
    ensemble='NVT',
    start_temp=start_temp,
    end_temp=end_temp,
    simulation_time = simulation_time,
    time_step = time_step,
    copy_vasp_outputs = True,
    supercell_matrix_kwargs: Optional[dict] = None
):
    """
    This workflow will perform AIMD on VASP,

    Args:
        structure: Initial structure.
        common_settings: Common settings dict. Supports "VASP_CMD", "DB_FILE",
            and "user_incar_settings" keys.
        vasp_input_set: Vasp input set for perturbed structure calculations.
        copy_vasp_outputs: Whether or not to copy previous vasp outputs.
        supercell_matrix_kwargs: Options that control the size of the supercell.
            Will be passed directly to CubicSupercellTransformation in
            pymatgen.transformations.advanced_transformations. Note, a diagonal
            supercell is required to calculate lattice thermal conductivity.
        ensemble: NVT, NPT, etc.
   """

    parent_structure = structure.as_dict()
    supercell_matrix_kwargs = supercell_matrix_kwargs or {}
    common_settings = _get_common_settings(common_settings)
    db_file = common_settings["DB_FILE"]

    logger.debug('Transforming supercell!')
    logger.debug('KWARGS: \n {}'.format(supercell_matrix_kwargs))
    st = CubicSupercellTransformation(**supercell_matrix_kwargs)
    supercell = st.apply_transformation(structure)
    supercell_matrix = st.transformation_matrix
    logger.debug('SUPERCELL MATRIX: \n {}'.format(supercell_matrix))

    nsteps_total = ceil(simulation_time/time_step)
    nsteps_each = 500
    nfws = ceil(nsteps_total/nsteps_each)
    parents=None

    vasp_input_set = ARCMDSet(
        structure=supercell,
        ensemble=ensemble,
        start_temp=start_temp,
        end_temp=end_temp,
        nsteps=nsteps_each,
        time_step=time_step,
        reciprocal_density=1,
        small_gap_multiply=None
    )
    for ifw in range(nfws):
        if ifw>0:
            vasp_input_set.update()
            parents = aimd_wf.fws[-1]
            start_temp = end_temp

        aimd_fw = ARCMDFW( start_structure,
                        start_temp,
                        end_temp,
                        nsteps_fw,
                        name="AIMD_segment_{}".format(ifw+1),
                        vasp_input_set=vasp_input_set,
                        vasp_cmd=VASP_CMD,
                        wall_time=1080000,
                        db_file=DB_FILE,
                        parents=parents,
                        copy_vasp_outputs=copy_vasp_outputs,
                        override_default_vasp_params=None,
                        **kwargs,)
        pass_task = pass_vasp_result(
            pass_dict={
            "parent_structure": parent_structure,
            "supercell_matrix": supercell_matrix,
            "forces": ">>output.ionic_steps.:-1.forces",
            "stress": ">>output.ionic_steps.:-1.stress",
            "configurations": ">>output.ionic_steps.:-1.structure",
            "final_configuration": ">>output.ionic_steps.-1.structure"
            },
            mod_spec_cmd="_push",
            mod_spec_key="aimd",
        )
        aimd_fw.tasks.append(pass_task)
        if ifw==0:
            wf = Workflow.from_Firework(aimd_fw)
        else:
            wf.append_wf(Workflow.from_Firework(aimd_fw), wf.fws[-1].fw_id)

    collect_fw = CollectMDSegmentsFW(
        db_file=db_file,
        parents=wf.fws[-nfws:]
    )
    wf.append_wf(
            Workflow.from_Firework(collect_fw), [fw.fw_id for fw in wf.fws[-nfws:]]
        )

    # Add workflow metadata to all relevant *ToDb tasks.
    wf_meta = {"wf_name": "AIMDWF", "wf_version": _WF_VERSION}
    for task_name in ["VaspToDb"]:
        wf = add_additional_fields_to_taskdocs(
            wf, {"wf_meta": wf_meta}, task_name_constraint=task_name
        )

    return wf


def _get_common_settings(common_settings: Optional[Dict]):
    common_settings = common_settings or {}
    for k, v in _DEFAULT_SETTINGS.items():
        if k not in common_settings:
            common_settings[k] = v

    user_incar_settings = deepcopy(_aimd_user_incar_settings)
    user_incar_settings.update(common_settings.get("user_incar_settings", {}))
    common_settings["user_incar_settings"] = user_incar_settings

    return common_settings


