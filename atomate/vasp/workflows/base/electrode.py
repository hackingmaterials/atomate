from typing import List

from fireworks import Workflow
from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher

from atomate.vasp.config import DB_FILE
from atomate.vasp.firetasks.electrode_tasks import AnalyzeChgcar, GetInsertionCalcs

__author__ = "Jimmy Shen"
__email__ = "jmmshn@lbl.gov"

from atomate.vasp.fireworks import Firework, OptimizeFW, StaticFW, pass_vasp_result
from atomate.common.powerups import powerup_by_kwargs

"""
Define workflow related to battery material simulation --- they all have a working ion
"""


def get_ion_insertion_wf(
    structure: Structure,
    working_ion: str,
    structure_matcher: StructureMatcher = None,
    db_file: str = DB_FILE,
    vasptodb_kwargs: dict = None,
    volumetric_data_type: str = "CHGCAR",
    vasp_powerups: List[dict] = None,
    max_insertions: int = 4,
    allow_fizzled_parents: bool = True,
    optimizefw_kwargs: dict = None,
    staticfw_kwargs: dict = None,
):
    """
    Take the output static worflow and iteratively insert working ions based on charge density analysis.

    The workflow performs the following tasks.
    (StaticFW) <- Recieved dat inserted task_id from this workflow
    (AnalyzeChgcar) <- Obtain the set of possible unique insertions using the stored charge density
    (GetInsertionCalcs) <- This task contains the dynamic workflow creation that will keep inserting working ions

    Args:
        structure: The host structure to begin inserting on
        working_ion: The working ion to be inserted at each step
        structure_matcher: StructureMatcher object used to define topotactic insertion
        db_file: The db_file that defines the VASP output database
        vasptodb_kwargs: vasptodb_kwargs for the static workflow
        volumetric_data_type: the type of volumetric data used to determine the insertion sites
        vasp_powerups: additional powerups given to all the dynamically created workflows
        optimizefw_kwargs: additional kwargs for all the OptimizeFWs
        staticfw_kwargs: additional kwargs for all the StaticFWs
    """

    if not structure.is_ordered:
        raise ValueError(
            "Please obtain an ordered approximation of the input structure."
        )

    if optimizefw_kwargs is None:
        optimizefw_kwargs = {}
    if staticfw_kwargs is None:
        staticfw_kwargs = {}

    # Configured the optimization and static FWs for the base material
    vasptodb_kwargs = vasptodb_kwargs if vasptodb_kwargs is not None else {}
    vasptodb_kwargs_vol_data = {"CHGCAR": ["CHGCAR"], "AECCAR": ["AECCAR0", "AECCAR2"]}

    vasptodb_kwargs.update(
        {
            "store_volumetric_data": vasptodb_kwargs_vol_data[volumetric_data_type],
            "task_fields_to_push": {"base_task_id": "task_id"},
        }
    )

    opt_wf = OptimizeFW(structure=structure, db_file=db_file, **optimizefw_kwargs)

    pass_task = pass_vasp_result(
        filename="vasprun.xml.relax2.gz",
        pass_dict=">>output.ionic_steps.-1.structure",
        mod_spec_key="prev_calc_structure",
    )
    opt_wf.tasks.append(pass_task)

    static_wf = StaticFW(
        structure=structure,
        vasptodb_kwargs=vasptodb_kwargs,
        db_file=db_file,
        parents=[opt_wf],
        spec_structure_key="prev_calc_structure",
        **staticfw_kwargs
    )

    wf_name = "{}-{}".format(
        structure.composition.reduced_formula if structure else "unknown",
        "insertion",
    )

    # Configure the analysis FW
    analysis_wf = Firework(
        [AnalyzeChgcar(), GetInsertionCalcs()],
        parents=[static_wf],
        name="Charge Density Analysis-0",
    )
    analysis_wf.spec["working_ion"] = working_ion

    # Crate the initial workflow
    wf = Workflow([opt_wf, static_wf, analysis_wf], name=wf_name)

    # Apply the vasp powerup if present
    if vasp_powerups is not None:
        wf = powerup_by_kwargs(wf, vasp_powerups)
        for fw in wf.fws:
            fw.spec["vasp_powerups"] = vasp_powerups

    if structure_matcher is not None:
        sm_dict = structure_matcher.as_dict()
        for fw in wf.fws:
            fw.spec["structure_matcher"] = sm_dict

    # write the persistent specs to all the fws
    # Note this is probably redundant but is easier
    for fw in wf.fws:
        fw.spec["db_file"] = db_file
        fw.spec["max_insertions"] = max_insertions
        fw.spec["vasptodb_kwargs"] = vasptodb_kwargs
        fw.spec["staticfw_kwargs"] = staticfw_kwargs
        fw.spec["optimizefw_kwargs"] = optimizefw_kwargs
        fw.spec["allow_fizzled_parents"] = allow_fizzled_parents
        fw.spec["volumetric_data_type"] = volumetric_data_type

    return wf
