# coding: utf-8


from pymatgen.io.vasp.sets import MPScanRelaxSet

from atomate.vasp.config import ADD_WF_METADATA

from atomate.vasp.powerups import (
    add_wf_metadata,
    add_common_powerups,
)
from atomate.vasp.workflows.base.core import get_wf

__author__ = "Ryan Kingsbury, Shyam Dwaraknath, Anubhav Jain"
__email__ = "rkingsbury@lbl.gov, shyamd@lbl.gov, ajain@lbl.gov"


def wf_scan_opt(structure, c=None):
    """
    Structure optimization using the SCAN metaGGA functional.

    This workflow performs a 2-step optmization. The first step
    is a conventional GGA run and serves to precondition the geometry and
    wavefunctions. The second step is a SCAN structure optimization.

    The first optimization is force converged with EDIFFG = -0.05,
    and the second optimization is force converged with EDIFFG=-0.02.
    """

    c = c or {}
    user_incar_settings = c.get("USER_INCAR_SETTINGS", {})
    vdw = c.get("vdw")
    bandgap = c.get("bandgap", 0)

    wf = get_wf(
        structure,
        "SCAN_optimization.yaml",
        vis=MPScanRelaxSet(
            structure,
            user_incar_settings=user_incar_settings,
            vdw=vdw,
            bandgap=bandgap
        ),
        params=[{"name": "SCAN optimization"}],
    )

    wf = add_common_powerups(wf, c)

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

    return wf
