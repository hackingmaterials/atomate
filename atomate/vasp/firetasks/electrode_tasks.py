import math
from collections import defaultdict

from fireworks import FiretaskBase, explicit_serialize, FWAction, Firework, Workflow
from pymatgen import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher

from atomate.utils.utils import env_chk, get_logger
from pymatgen.analysis.defects.utils import ChargeInsertionAnalyzer
from atomate.vasp.fireworks.core import OptimizeFW, StaticFW

__author__ = "Jimmy Shen"
__email__ = "jmmshn@lbl.gov"

from atomate.vasp.database import VaspCalcDb
from atomate.vasp.firetasks import pass_vasp_result
from atomate.vasp.powerups import powerup_by_kwargs, POWERUP_NAMES

logger = get_logger(__name__)

sm = StructureMatcher()


@explicit_serialize
class AnalyzeChgcar(FiretaskBase):
    """
    Process the charge density from a particular task.
    fw_spec requires:
    - "base_task_id" : the task_id for the charge density that is being analyzed
    - "ChargeInsertionAnalyzer_kwargs": The kwargs to overwrite default values in the ChargeInsertionAnalyzer
    - "db_file": The db_file that instantiates the tasks database
    """

    required_params = ["base_task_id"]

    def run_task(self, fw_spec):
        base_task_id = self.get("base_task_id")
        cia_kwargs = fw_spec.get("ChargeInsertionAnalyzer_kwargs", dict())

        # get the database connection
        db_file = env_chk(self.get("db_file"), fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        chgcar = mmdb.get_chgcar(task_id=base_task_id)

        cia = ChargeInsertionAnalyzer(chgcar, **cia_kwargs)
        cia.get_labels()

        insert_sites = []

        for itr, li_site in cia._extrema_df.iterrows():
            li_site = self._extrema_df.iloc[itr]
            insert_sites.append([li_site["a"], li_site["b"], li_site["c"]])

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
class GetNewCalcs(FiretaskBase):
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
        base_structure = fw_spec.get("base_structure", None)
        working_ion = fw_spec.get("working_ion", "Li")

        if base_structure is None:
            raise RuntimeError(
                "No base structure was passed to generate new insertion calculations."
            )

        base_structure = Structure.from_dict(base_structure)
        new_fws = []

        for isite in insert_sites:
            inserted_structure = base_structure.copy()
            fpos = [isite["a"], isite["b"], isite["c"]]
            inserted_structure.insert(0, working_ion, fpos, properties={"magmom": 0})

            additional_fields = {"insertion_fpos": fpos, "base_task_id": base_task_id}

            # Create new fw
            wf = Workflow(OptimizeFW(inserted_structure))
            wf = get_powereup_wf(wf, fw_spec, additional_fields=additional_fields)

            pass_dict = {
                "structure": ">>output.ionic_steps.-1.structure",
            }

            pass_task = pass_vasp_result(
                pass_dict=pass_dict,
                mod_spec_cmd="_push",
                mod_spec_key="inserted_tasks",
            )

            wf.tasks.append(pass_task)
            new_fws.append(wf)

        if len(new_fws) == 0:
            return

        # analysis workflow that returns a
        check_task = (
            CollectInsertedCalcs()
        )  # Ask Alex: Does this just grab the "host_structure" key
        check_fw = Firework([check_task], parents=new_fws)  # Allow fizzled parent
        check_fw.spec["_allow_fizzled_parents"] = True

        return FWAction(additions=new_fws + [check_fw])


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
        reference_struct = Structure.from_dict(fw_spec.get("host_structure"))
        results = fw_spec.get("inserted_results", [])
        working_ion = fw_spec.get("working_ion", "Li")
        sm_dict = fw_spec.get("StructureMatcher", {})

        sm_dict["ignore_species"] = [working_ion]

        try:
            sm = StructureMatcher.from_dict(sm_dict)
        except KeyError:
            sm = Structure(ignore_species=[working_ion])

        if len(results) == 0:
            # can happen if all parents fizzled
            raise RuntimeError("No insertion calculation completed")

        best_res = {"energy": math.inf}

        n_completed = 0
        for ires in results:
            n_completed += 1
            ires_struct_ = Structure.from_dict(ires["structure"])
            if (
                sm.fit(ires_struct_, reference_struct)
                and ires["energy"] < best_res["energy"]
            ):
                best_res = ires

        if "task_id" not in best_res:
            # No matching structure was found in the completed results
            if n_completed > 0:
                return FWAction()  # TODO maybe deffuse children here?

        # Get the new structures
        return FWAction(update_spec=[{"base_task_id": best_res["task_id"]}])


@explicit_serialize
class SubmitMostStable(FiretaskBase):
    """
    For the best structure submit the Static WF
    """

    _fw_name = "SubmitBestInsertion"

    def run_task(self, fw_spec):
        inserted_structure = fw_spec.get("inserted_structure")
        inserted_structure = Structure.from_dict(inserted_structure)

        wf = Workflow(
            StaticFW(structure=inserted_structure)
        )  # how to set the structure here
        wf.tasks.append(AnalyzeChgcar())
        wf.tasks.append(GetNewCalcs())
        wf = get_powereup_wf(wf, fw_spec)
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
        d_pu["add_additional_fields_to_taskdocs"].update(additional_fields)
    p_kwargs = {k: d_pu[k] for k in POWERUP_NAMES if k in d_pu}
    return powerup_by_kwargs(wf, **p_kwargs)
