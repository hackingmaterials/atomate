import math
from collections import defaultdict

from fireworks import FiretaskBase, explicit_serialize, FWAction, Firework, Workflow
from pymatgen import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher

from atomate.utils.utils import env_chk, get_logger
from pymatgen.analysis.defects.utils import ChargeInsertionAnalyzer

from atomate.vasp.config import DB_FILE
from atomate.vasp.fireworks.core import OptimizeFW, StaticFW

from atomate.vasp.database import VaspCalcDb
from atomate.vasp.firetasks import pass_vasp_result
from atomate.vasp.powerups import powerup_by_kwargs, POWERUP_NAMES

__author__ = "Jimmy Shen"
__email__ = "jmmshn@lbl.gov"

logger = get_logger(__name__)

sm = StructureMatcher()


# Helper function to allow specific keys to be propagated to newly created workflows
def update_wf_keys(wf, fw_spec):
    PASS_KEYS = [
        "working_ion",
        "db_file",
        "vasp_powerups",
        "base_task_id",
        "max_insertions",
        "base_structure",
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
    - "db_file": The db_file that instantiates the tasks database
    - "max_insertions": Restrict the number of insertion sites based on charge density (default 5)
    """

    optional_params = ["db_file"]

    def run_task(self, fw_spec):
        max_insertions = fw_spec.get("max_insertions", 5)
        base_task_id = fw_spec.get("base_task_id")
        logger.info(f"Identifying sites for task : {base_task_id}")

        cia_kwargs = fw_spec.get("ChargeInsertionAnalyzer_kwargs", dict())

        # get the database connection
        fw_spec["db_file"] = DB_FILE
        db_file = env_chk(fw_spec.get("db_file"), fw_spec)
        logger.info(f"DB_FILE: {db_file}")
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        chgcar = mmdb.get_chgcar(task_id=base_task_id)

        cia = ChargeInsertionAnalyzer(chgcar, **cia_kwargs)
        cia.get_labels()

        insert_sites = []
        seent = set()

        cia._extrema_df.sort_values(by=["avg_charge_den"], inplace=True)
        for itr, li_site in cia._extrema_df.iterrows():
            if len(insert_sites) >= max_insertions:
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
    """

    def run_task(self, fw_spec):

        insert_sites = fw_spec.get("insert_sites", [])
        base_task_id = fw_spec.get("base_task_id")
        base_structure = Structure.from_dict(fw_spec.get("base_structure"))
        working_ion = fw_spec.get("working_ion")

        new_fws = []
        for itr, isite in enumerate(insert_sites):
            inserted_structure = base_structure.copy()
            fpos = isite
            inserted_structure.insert(0, working_ion, fpos, properties={"magmom": None})

            additional_fields = {"insertion_fpos": fpos, "base_task_id": base_task_id}

            # Create new fw
            fw = OptimizeFW(inserted_structure, name=f"insert-opt-{itr}")
            fw.tasks[-1]["additional_fields"].update(additional_fields)

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
            name="Collect Inserted Calcs",
        )  # Allow fizzled parent
        check_fw.spec["_allow_fizzled_parents"] = True

        wf = Workflow(new_fws + [check_fw])
        wf = get_powereup_wf(wf, fw_spec)

        update_wf_keys(wf, fw_spec)
        return FWAction(additions=[wf])


@explicit_serialize
class CollectInsertedCalcs(FiretaskBase):
    """
    Aggregate the outputs of the inserted calculations, take the lowest energy completed calculation
    and pass it on for charge density analysis

    Required keys in fw_spec:
    - inserted_results: list of parsed results from relaxation calcs each result needs to be:
            {task_id : id, energy: -123.45, structure: struct}
    - StructureMatcher: as_dict of a pymatgen structure matcher
    - working_ion: name of the working ion
    """

    _fw_name = "CollectInsertedCalc"

    def run_task(self, fw_spec):
        reference_struct = fw_spec.get("base_structure")
        results = fw_spec.get("inserted_results", [])
        working_ion = fw_spec.get("working_ion")
        sm_dict = fw_spec.get("StructureMatcher", {})

        sm_dict["ignore_species"] = [working_ion]

        try:
            sm = StructureMatcher.from_dict(sm_dict)
        except KeyError:
            sm = StructureMatcher(ignored_species=[working_ion])

        if len(results) == 0:
            # can happen if all parents fizzled
            raise RuntimeError("No insertion calculation completed")

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
                return FWAction(
                    defuse_children=True
                )  # TODO maybe diffuse children here?

        # Get the new structures
        return FWAction(update_spec={"optimal_structure": best_res["structure"]})


@explicit_serialize
class SubmitMostStable(FiretaskBase):
    """
    For the best structure submit the Static WF
    """

    _fw_name = "SubmitBestInsertion"

    def run_task(self, fw_spec):
        inserted_structure = fw_spec.get("optimal_structure")

        vasptodb_kwargs = {
            "store_volumetric_data": ["CHGCAR"],
            "task_fields_to_push": {"base_task_id": "task_id"},
        }

        fw1 = StaticFW(
            inserted_structure, vasptodb_kwargs=vasptodb_kwargs, db_file=DB_FILE
        )
        fw2 = Firework(
            [AnalyzeChgcar(), GetInsertionCalcs()],
            name="Charge Density Analysis",
            parents=fw1,
        )
        wf = Workflow([fw1, fw2], name="Obtain inserted CHG")
        wf = get_powereup_wf(wf, fw_spec)
        update_wf_keys(wf, fw_spec)
        return FWAction(additions=[wf])


def get_powereup_wf(wf, fw_spec, additional_fields=None):
    """
    Check the fw_spec['vasp_powerups'] for powerups and apply them to a workflow.
    Add/overwrite the additional fields in the fw_spec with user inputs
    Args:
        fw_spec: the fw_spec of the current workflow
        additional_fields: The additional fields to be added to the task document
    Returns:
        Updated workflow
    """
    d_pu = defaultdict(dict)
    d_pu.update(fw_spec.get("vasp_powerups", {}))
    if additional_fields is not None:
        d_pu["add_additional_fields_to_taskdocs"].update(
            {"update_dict": additional_fields}
        )
    p_kwargs = {k: d_pu[k] for k in POWERUP_NAMES if k in d_pu}
    return powerup_by_kwargs(wf, **p_kwargs)
