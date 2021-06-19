import math

from fireworks import FiretaskBase, explicit_serialize, FWAction, Firework, Workflow
from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher

from atomate.utils.utils import env_chk, get_logger
from pymatgen.analysis.defects.utils import ChargeInsertionAnalyzer

from atomate.vasp.config import DB_FILE
from atomate.vasp.fireworks.core import OptimizeFW, StaticFW

from atomate.vasp.database import VaspCalcDb
from atomate.vasp.firetasks import pass_vasp_result
from atomate.common.powerups import powerup_by_kwargs

__author__ = "Jimmy Shen"
__email__ = "jmmshn@lbl.gov"

"""
A recursive list of Tasks to keep inserting ions into a structure.
Should only be used by the `get_ion_insertion_wf` function.

Note:
The workflow passes data to fw_spec extensively and requires specific fields in the spec to be updated.
Example, the base_taks_id must be stored and the spec and updated as the workflow runs so you have to set
```
{
    "store_volumetric_data": vasptodb_kwargs_vol_data[volumetric_data_type],
    "task_fields_to_push": {"base_task_id": "task_id"},
}
```
In your vasptodb_kwargs.  Which is currently done in `get_ion_insertion_wf`
"""

logger = get_logger(__name__)

SM_DICT = StructureMatcher().as_dict()


def update_wf_keys(wf, fw_spec):
    """
    Helper function to allow specific keys to be propagated to newly created workflows
    update the wf in place with the keys in the fw_spec dictionary
    Args:
        wf: workflow to modify
        fw_spec: dictionary to be passed
    """

    PASS_KEYS = [
        "working_ion",
        "db_file",
        "vasp_powerups",
        "base_task_id",
        "attempt_insertions",
        "max_inserted_atoms",
        "base_structure",
        "vasptodb_kwargs",
        "staticfw_kwargs",
        "optimizefw_kwargs",
        "structure_matcher",
        "allow_fizzled_parents",
        "volumetric_data_type",
    ]
    for k in PASS_KEYS:
        if k in fw_spec:
            for fw in wf.fws:
                fw.spec.update({k: fw_spec.get(k, None)})


@explicit_serialize
class AnalyzeChgcar(FiretaskBase):
    """
    Process the charge density from a particular task.
    fw_spec requires:
    - "base_task_id" : the task_id for the charge density that is being analyzed
    - "ChargeInsertionAnalyzer_kwargs": The kwargs to overwrite default values in the ChargeInsertionAnalyzer
    - "volumetric_data_type": the type of Charge density file to be used for analysis "CHGCAR"/"AECCAR"
    fw_spec optional:
    - "attempt_insertions": Restrict the number of insertion sites based on charge density (default 4)
    """

    def run_task(self, fw_spec):
        attempt_insertions = fw_spec.get("attempt_insertions", 4)
        volumetric_data_type = fw_spec.get("volumetric_data_type")
        base_task_id = fw_spec.get("base_task_id")

        logger.info(f"Identifying sites for task : {base_task_id}")

        cia_kwargs = fw_spec.get("ChargeInsertionAnalyzer_kwargs", dict())

        # get the database connection
        db_file = env_chk(DB_FILE, fw_spec)
        logger.info(f"DB_FILE: {db_file}")
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        if volumetric_data_type == "CHGCAR":
            chgcar = mmdb.get_chgcar(task_id=base_task_id)
        elif volumetric_data_type == "AECCAR":
            chgcar = mmdb.get_aeccar(task_id=base_task_id)
            chgcar = chgcar["aeccar0"] + chgcar["aeccar2"]

        cia = ChargeInsertionAnalyzer(chgcar, **cia_kwargs)
        cia.get_labels()

        insert_sites = []
        seent = set()

        cia._extrema_df.sort_values(by=["avg_charge_den"], inplace=True)
        for itr, li_site in cia._extrema_df.iterrows():
            if len(insert_sites) >= attempt_insertions:
                break
            li_site = cia._extrema_df.iloc[itr]
            lab = li_site["site_label"]
            if lab not in seent:
                insert_sites.append([li_site["a"], li_site["b"], li_site["c"]])
                seent.add(lab)

        logger.info(
            f"Found {len(insert_sites)} insertion sites for task : {base_task_id}"
        )

        # Since we are only analyzing a single charge density
        # The we will only perform one set of atoms insertions at a given time
        # Thus, we just have to update the global fw_spec
        # (no need to pass this information from fw to fw)
        return FWAction(
            update_spec={
                "insert_sites": insert_sites,
                "base_task_id": base_task_id,
                "base_structure": chgcar.structure.as_dict(),
            }
        )


@explicit_serialize
class GetInsertionCalcs(FiretaskBase):
    """
    Using a list of insertion sites and a base structure, get a new set of calculations

    fw_spec requires:
    - "base_task_id" : the task_id for the charge density that is being analyzed
    - "base_structure": the base structure that the ion will be inserted on
    - "insert_sites": The list of insertion sites
    - "working_ion": The name of the ion to be inserted

    Optional fw_spec:
    - "allow_fizzled_parents": If True, perform the static calculation and charge density
        analysis as long as a single optimization is finished. Default is False
    - "optimizefw_kwargs":  kwargs for the OptimizeFW created by the present task
    """

    def run_task(self, fw_spec):

        insert_sites = fw_spec.get("insert_sites", [])
        base_task_id = fw_spec.get("base_task_id")
        base_structure = Structure.from_dict(fw_spec.get("base_structure"))
        working_ion = fw_spec.get("working_ion")
        allow_fizzled_parents = fw_spec.get("allow_fizzled_parents", False)
        optimizefw_kwargs = fw_spec.get("optimizefw_kwargs", {})
        n_ion = int(base_structure.composition.element_composition[working_ion]) + 1

        new_fws = []
        for itr, isite in enumerate(insert_sites):
            inserted_structure = base_structure.copy()
            fpos = isite
            inserted_structure.insert(0, working_ion, fpos, properties={"magmom": 1.0})
            additional_fields = {"insertion_fpos": fpos, "base_task_id": base_task_id}

            # Create new fw
            fw = OptimizeFW(
                inserted_structure,
                name=f"structure optimization-{itr}",  # structure are rough guesses
                **optimizefw_kwargs,
            )
            fw.tasks[-1]["additional_fields"].update(
                additional_fields
            )  # easy way to update a field

            pass_dict = {
                "structure": ">>output.ionic_steps.-1.structure",
                "energy": "a>>final_energy",
            }

            pass_task = pass_vasp_result(
                filename="vasprun.xml.relax2.gz",
                pass_dict=pass_dict,
                mod_spec_cmd="_push",
                mod_spec_key="inserted_results",
            )

            fw.tasks.append(pass_task)
            new_fws.append(fw)

        if len(new_fws) == 0:
            return

        check_fw = Firework(
            [CollectInsertedCalcs(), SubmitMostStable()],
            parents=new_fws,
            name=f"Collect Inserted Calcs-{n_ion}",
        )
        # Allow fizzled parent
        check_fw.spec["_allow_fizzled_parents"] = allow_fizzled_parents

        wf = Workflow(new_fws + [check_fw])
        wf = get_powerup_wf(wf, fw_spec)

        update_wf_keys(wf, fw_spec)
        return FWAction(additions=[wf])


@explicit_serialize
class CollectInsertedCalcs(FiretaskBase):
    """
    Aggregate the outputs of the inserted calculations, take the lowest energy completed calculation
    and pass it on for charge density analysis

    Required keys in fw_spec:
    - "base_structure": the base structure that the ion will be inserted on
    - "inserted_results": list of parsed results from relaxation calcs each result needs to be:
            {task_id : id, energy: -123.45, structure: struct}
    - "working_ion": name of the working ion

    Optional keys in fw_spec:
    - "structure_matcher": as_dict of a pymatgen structure matcher
    - "max_inserted_atoms": the maximum allowed number of inserted atoms
    - "staticfw_kwargs": as_dict of a pymatgen structure matcher

    """

    _fw_name = "CollectInsertedCalc"

    def run_task(self, fw_spec):
        reference_struct = fw_spec.get("base_structure")
        results = fw_spec.get("inserted_results", [])
        max_inserted_atoms = fw_spec.get("max_inserted_atoms", None)
        working_ion = fw_spec.get("working_ion")

        sm_dict = fw_spec.get("structure_matcher", SM_DICT)
        if not isinstance(sm_dict, dict):
            sm_dict = sm_dict.as_dict()
        sm_dict["ignored_species"] = [working_ion]
        sm = StructureMatcher.from_dict(sm_dict)

        if len(results) == 0:
            # can happen if all parents fizzled
            raise RuntimeError("No insertion calculation completed")

        tmp_struct_ = results[0]["structure"]  # type: Structure
        n_ions = int(tmp_struct_.composition.as_dict()[working_ion])
        if max_inserted_atoms is not None and n_ions >= max_inserted_atoms:
            return FWAction(defuse_children=True)

        best_res = {"energy": math.inf}

        n_completed = 0
        for ires in results:
            n_completed += 1
            ires_struct_ = ires["structure"]
            if (
                sm.fit(ires_struct_, reference_struct)
                and ires["energy"] < best_res["energy"]
            ):
                logger.info("Found new optimal_structure")
                best_res = ires

        if "structure" not in best_res:
            # No matching structure was found in the completed results
            if n_completed > 0:
                return FWAction(defuse_children=True)

        # Get the new structures
        return FWAction(update_spec={"optimal_structure": best_res["structure"]})


@explicit_serialize
class SubmitMostStable(FiretaskBase):
    """
    For the best structure submit the Static WF

    Required keys in fw_spec:
    - optimal_structure: the most stable relaxed structure
    - working_ion: name of the working ion
    - staticfw_kwargs = fw_spec.get("staticfw_kwargs")
    """

    _fw_name = "SubmitBestInsertion"

    def run_task(self, fw_spec):
        inserted_structure = fw_spec.get("optimal_structure")
        working_ion = fw_spec.get("working_ion")
        vasptodb_kwargs = fw_spec.get("vasptodb_kwargs")
        staticfw_kwargs = fw_spec.get("staticfw_kwargs", {})

        fw1 = StaticFW(
            inserted_structure,
            vasptodb_kwargs=vasptodb_kwargs,
            db_file=DB_FILE,
            **staticfw_kwargs,
        )
        n_ion = int(inserted_structure.composition.element_composition[working_ion])
        fw2 = Firework(
            [AnalyzeChgcar(), GetInsertionCalcs()],
            name=f"Charge Density Analysis-{n_ion}",
            parents=fw1,
        )
        wf = Workflow([fw1, fw2], name=f"Obtain inserted sites-{n_ion}")
        wf = get_powerup_wf(wf, fw_spec)
        update_wf_keys(wf, fw_spec)
        return FWAction(additions=[wf])


def get_powerup_wf(wf, fw_spec, additional_fields=None):
    """
    Check the fw_spec['vasp_powerups'] for powerups and apply them to a workflow.
    Add/overwrite the additional fields in the fw_spec with user inputs
    Args:
        fw_spec: the fw_spec of the current workflow
        additional_fields: The additional fields to be added to the task document
    Returns:
        Updated workflow
    """
    powerup_list = []
    powerup_list.extend(fw_spec.get("vasp_powerups", []))
    if additional_fields is not None:
        powerup_list.append(
            {
                "powerup_name": "add_additional_fields_to_taskdocs",
                "kwargs": {"update_dict": additional_fields},
            }
        )
    return powerup_by_kwargs(wf, powerup_list)
