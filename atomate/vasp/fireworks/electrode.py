from atomate.vasp.firetasks.electrode_tasks import AnalyzeChgcar

__author__ = "Jimmy Shen"
__email__ = "jmmshn@lbl.gov"

from atomate.vasp.fireworks import StaticFW

"""
Define worflows related to battery material simulation --- they all have a working ion
"""


class IonInsertionWF(StaticFW):
    """
    Take the output static worflow and iteratively insert working ions based on charge density analysis.

    The workflow performs the following tasks.
    (StaticFW) <- Recieved dat inserted task_id from this workflow
    (AnalyzeChgcar)
    (GetInsertionCalcs) <- This task contains the dynamic workflow creation that will keep inserting working ions
    """

    def __init__(
        self,
        structure,
        vasptodb_kwargs: dict = None,
        volumetric_data_type: str = "CHGCAR",
        vasp_powerups: dict = None,
        **kwargs
    ):
        # set up the StaticFW and make sure the base_task_id field is populated
        vasptodb_kwargs = vasptodb_kwargs if vasptodb_kwargs is not None else {}
        vasptodb_kwargs.update(
            {
                "store_volumetric_data": [volumetric_data_type],
                "task_fields_to_push": {"base_task_id": "task_id"},
            }
        )

        super().__init__(structure=structure, vasptodb_kwargs=vasptodb_kwargs)

        if vasp_powerups is not None:
            for fw in self.fws:
                fw.spec["vasp_powerups"] = vasp_powerups

        self.tasks.extend([AnalyzeChgcar()])
