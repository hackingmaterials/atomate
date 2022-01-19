"""
Define workflow related to battery material simulation --- they all have a working ion
"""
from typing import List

from fireworks import Firework, Workflow
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure

from atomate.common.powerups import powerup_by_kwargs
from atomate.vasp.config import DB_FILE
from atomate.vasp.firetasks import pass_vasp_result
from atomate.vasp.firetasks.electrode_tasks import AnalyzeChgcar, GetInsertionCalcs
from atomate.vasp.fireworks import OptimizeFW, StaticFW

__author__ = "Jimmy Shen"
__email__ = "jmmshn@lbl.gov"


def get_ion_insertion_wf(
    structure: Structure,
    working_ion: str,
    structure_matcher: StructureMatcher = None,
    db_file: str = DB_FILE,
    vasptodb_kwargs: dict = None,
    volumetric_data_type: str = "CHGCAR",
    vasp_powerups: List[dict] = None,
    attempt_insertions: int = 4,
    max_inserted_atoms: int = None,
    allow_fizzled_parents: bool = True,
    optimizefw_kwargs: dict = None,
    staticfw_kwargs: dict = None,
):
    """
    Take the output static workflow and iteratively insert working ions based on charge density analysis.

    The workflow performs the following tasks.
    (StaticFW) <- Received dat inserted task_id from this workflow
    (AnalyzeChgcar) <- Obtain the set of possible unique insertions using the stored charge density
    (GetInsertionCalcs) <- This task contains the dynamic workflow creation that will keep inserting working ions

    If you use this workflow please cite the following paper:
        Shen, J.-X., Horton, M., & Persson, K. A. (2020).
        A charge-density-based general cation insertion algorithm for generating new Li-ion cathode materials.
        npj Comput. Mater., 6(161), 1–7. doi: 10.1038/s41524-020-00422-3

    Args:
        structure: The host structure to begin inserting on
        working_ion: The working ion to be inserted at each step
        structure_matcher: StructureMatcher object used to define topotactic insertion
        db_file: The db_file that defines the VASP output database
        vasptodb_kwargs: vasptodb_kwargs for the static workflow
        volumetric_data_type: the type of volumetric data used to determine the insertion sites
        vasp_powerups: additional powerups given to all the dynamically created workflows
        optimizefw_kwargs: additional kwargs for all the OptimizeFWs
        max_inserted_atoms: the limit on the total number of ions to insert
        attempt_insertions: number of insertions sites to run at each insertion step
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

    # Create the initial workflow
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
        fw.spec["attempt_insertions"] = attempt_insertions
        fw.spec["vasptodb_kwargs"] = vasptodb_kwargs
        fw.spec["staticfw_kwargs"] = staticfw_kwargs
        fw.spec["optimizefw_kwargs"] = optimizefw_kwargs
        fw.spec["allow_fizzled_parents"] = allow_fizzled_parents
        fw.spec["volumetric_data_type"] = volumetric_data_type
        fw.spec["max_inserted_atoms"] = max_inserted_atoms

    return wf
