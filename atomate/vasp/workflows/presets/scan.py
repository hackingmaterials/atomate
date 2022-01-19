from atomate.vasp.config import ADD_WF_METADATA
from atomate.vasp.powerups import add_common_powerups, add_wf_metadata
from atomate.vasp.workflows.base.core import get_wf

__author__ = "Ryan Kingsbury, Shyam Dwaraknath, Anubhav Jain"
__email__ = "rkingsbury@lbl.gov, shyamd@lbl.gov, ajain@lbl.gov"


def wf_r2scan_opt(structure, c=None):
    """
    Structure optimization using the R2SCAN metaGGA functional.

    This workflow performs a 2-step optimization. The first step
    is a GGA structure optimization using the PBESol functional that serves to
    precondition the geometry and charge density. The second step is a
    R2SCAN structure optimization.

    The first optimization is force converged with EDIFFG = -0.05,
    and the second optimization is force converged with EDIFFG=-0.02.

    The bandgap from the first step is used to update the KSPACING parameter,
    which sets the appropriate number of k-points for the subsequent R2SCAN
    calculation.
    """

    c = c or {}
    vasp_input_set_params = {}
    if c.get("USER_INCAR_SETTINGS"):
        vasp_input_set_params["user_incar_settings"] = c.get("USER_INCAR_SETTINGS")

    if c.get("vdw"):
        vasp_input_set_params["vdw"] = c.get("vdw")

    if c.get("bandgap"):
        vasp_input_set_params["bandgap"] = c.get("bandgap")

    wf = get_wf(
        structure,
        "metagga_optimization.yaml",
        common_params={"vasp_input_set_params": vasp_input_set_params},
    )

    wf = add_common_powerups(wf, c)

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

    return wf


def wf_scan_opt(structure, c=None):
    """
    Structure optimization using the SCAN metaGGA functional.

    This workflow performs a 2-step optimization. The first step
    is a GGA structure optimization using the PBESol functional that serves to
    precondition the geometry and charge density. The second step is a
    SCAN structure optimization.

    The first optimization is force converged with EDIFFG = -0.05,
    and the second optimization is force converged with EDIFFG=-0.02.

    The bandgap from the first step is used to update the KSPACING parameter,
    which sets the appropriate number of k-points for the subsequent SCAN
    calculation.
    """

    c = c or {}
    # override the default R2SCAN functional with SCAN
    vasp_input_set_params = {"user_incar_settings": {"METAGGA": "SCAN"}}
    if c.get("USER_INCAR_SETTINGS"):
        vasp_input_set_params["user_incar_settings"] = c.get("USER_INCAR_SETTINGS")
        # if the user has supplied METAGGA, respect that setting instead
        if not vasp_input_set_params["user_incar_settings"].get("METAGGA"):
            vasp_input_set_params["user_incar_settings"]["METAGGA"] = "SCAN"

    if c.get("vdw"):
        vasp_input_set_params["vdw"] = c.get("vdw")

    if c.get("bandgap"):
        vasp_input_set_params["bandgap"] = c.get("bandgap")

    wf = get_wf(
        structure,
        "metagga_optimization.yaml",
        # override the default FW name to reflect the SCAN functional
        params=[{}, {"name": "SCAN structure optimization"}],
        common_params={"vasp_input_set_params": vasp_input_set_params},
    )

    wf = add_common_powerups(wf, c)

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

    return wf
