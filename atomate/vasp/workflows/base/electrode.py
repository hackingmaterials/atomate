from collections import Mapping

from fireworks import Workflow
from pymatgen import Structure

from atomate.vasp.config import DB_FILE
from atomate.vasp.firetasks.electrode_tasks import AnalyzeChgcar, GetInsertionCalcs

__author__ = "Jimmy Shen"
__email__ = "jmmshn@lbl.gov"

from atomate.vasp.fireworks import StaticFW, Firework
from atomate.vasp.powerups import powerup_by_kwargs

"""
Define worflows related to battery material simulation --- they all have a working ion
"""


class IonInsertionWF(StaticFW):
    def __init__(
        self,
        structure,
        vasptodb_kwargs: dict = None,
        volumetric_data_type: str = "CHGCAR",
        vasp_powerups: dict = None,
        **kwargs
    ):
        super().__init__(structure=structure, vasptodb_kwargs=vasptodb_kwargs)

        if vasp_powerups is not None:
            for fw in self.fws:
                fw.spec["vasp_powerups"] = vasp_powerups

        self.tasks.extend([AnalyzeChgcar()])


def get_ion_insertion_wf(
    structure: Structure,
    working_ion: str,
    db_file: str = DB_FILE,
    vasptodb_kwargs: dict = None,
    volumetric_data_type: str = "CHGCAR",
    vasp_powerups: dict = None,
    **kwargs
):
    """
    Take the output static worflow and iteratively insert working ions based on charge density analysis.

    The workflow performs the following tasks.
    (StaticFW) <- Recieved dat inserted task_id from this workflow
    (AnalyzeChgcar)
    (GetInsertionCalcs) <- This task contains the dynamic workflow creation that will keep inserting working ions

    Args:
        structure:
        working_ion:
        db_file:
        vasptodb_kwargs:
        volumetric_data_type:
        vasp_powerups:
        **kwargs:

    """

    if not structure.is_ordered:
        raise ValueError(
            "Please obtain an ordered approximation of the input structure."
        )

    # set up the StaticFW and make sure the base_task_id field is populated
    vasptodb_kwargs = vasptodb_kwargs if vasptodb_kwargs is not None else {}
    vasptodb_kwargs.update(
        {
            "store_volumetric_data": [volumetric_data_type],
            "task_fields_to_push": {"base_task_id": "task_id"},
        }
    )

    static_wf = StaticFW(
        structure=structure, vasptodb_kwargs=vasptodb_kwargs, db_file=db_file, **kwargs
    )

    wf_name = "{}-{}".format(
        structure.composition.reduced_formula if structure else "unknown", "insertion"
    )

    analysis_wf = Firework(
        [AnalyzeChgcar(), GetInsertionCalcs()],
        parents=[static_wf],
        name="Charge Density Analysis",
    )
    analysis_wf.spec["working_ion"] = working_ion
    wf = Workflow([static_wf, analysis_wf], name=wf_name)

    for fw in wf.fws:
        fw.spec["db_file"] = db_file

    vp_dict = {
        "add_modify_incar": {
            "modify_incar_params": {"incar_update": {"ISMEAR": 0, "SIGMA": 0.02}}
        }
    }

    if vasp_powerups is not None:
        rec_update(vp_dict, vasp_powerups)

    for fw in wf.fws:
        fw.spec["vasp_powerups"] = vp_dict

    wf = powerup_by_kwargs(wf, **vp_dict)

    return wf


def rec_update(d, u):
    for k, v in u.items():
        if isinstance(v, Mapping):
            d[k] = rec_update(d.get(k, {}), v)
        else:
            d[k] = v
    return d
