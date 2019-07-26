from fireworks import FiretaskBase, FWAction, explicit_serialize

from atomate.utils.utils import get_logger, env_chk
from atomate.vasp.database import VaspCalcDb

from pymatgen import Structure
from pymatgen.analysis.magnetism.heisenberg import HeisenbergMapper

__author__ = "Nathan C. Frey"
__email__ = "ncfrey@lbl.gov"

logger = get_logger(__name__)


@explicit_serialize
class HeisenbergModelMapping(FiretaskBase):
    """
    Map structures and energies to a Heisenberg model and compute exchange parameters for a given NN cutoff.

    Args:
        db_file (str): path to file containing the database credentials.
        exchange_wf_uuid (int): Unique id for record keeping.
    """

    required_params = [
        "db_file",
        "exchange_wf_uuid",
        "parent_structure",
        "cutoff",
        "tol",
    ]
    optional_params = []

    def run_task(self, fw_spec):
        db_file = env_chk(self["db_file"], fw_spec)
        wf_uuid = self["exchange_wf_uuid"]

        # Get magnetic orderings collection from db
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        mmdb.collection = mmdb.db["magnetic_orderings"]

        formula = self["parent_structure"].formula
        formula_pretty = self["parent_structure"].composition.reduced_formula

        # Get documents
        docs = list(
            mmdb.collection.find(
                {"wf_meta.wf_uuid": wf_uuid},
                ["task_id", "structure", "energy_per_atom"],
            )
        )

        # Get structures and energy / unit cell
        structures = [Structure.from_dict(s) for s in docs["structure"]]
        n_atoms = len(structures[0])
        energies = [e * n_atoms for e in docs["energy_per_atom"]]

        # Map system to a Heisenberg Model
        hmapper = HeisenbergMapper(structures, energies, self["cutoff"], self["tol"])

        # Get MSONable Heisenberg Model
        hmodel = hmapper.get_heisenberg_model()
        name = "heisenberg_model_" + str(self["cutoff"])
        task_doc = {
            "wf_meta.wf_uuid": wf_uuid,
            "wf_meta.analysis_task_id": docs["task_id"],
            "formula_pretty": formula_pretty,
            "nn_cutoff": self["cutoff"],
            "nn_tol": self["tol"],
            "heisenberg_model": hmodel,
        }

        # Exchange collection
        mmdb.collection = mmdb.db["exchange"]
        mmdb.collection.insert(task_doc)


@explicit_serialize
class HeisenbergConvergence(FiretaskBase):
    """
    Quick check to see if the Heisenberg model has "converged" for any particular nearest neighbor cutoff value in the sense that |J_ij| < |E_0| for all i, j. 

    If not, it doesn't make sense to do Monte Carlo and the system is either 1) not well described by the Heisenberg Model (perhaps the magnetic moments are not localized) or 2) not correctly captured by the provided low-energy magnetic orderings.

    Args:
        db_file (str): path to file containing the database credentials.
        exchange_wf_uuid (int): Unique id for record keeping.

    """

    required_params = ["db_file", "exchange_wf_uuid", "parent_structure"]
    optional_params = []

    def run_task(self, fw_spec):

        db_file = env_chk(self["db_file"], fw_spec)
        wf_uuid = self["exchange_wf_uuid"]

        # Get Heisenberg models from db
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        mmdb.collection = mmdb.db["exchange"]

        formula_pretty = self["parent_structure"].composition.reduced_formula

        # Get documents
        docs = list(
            mmdb.collection.find(
                {"wf_meta.wf_uuid": wf_uuid},
                ["task_id", "heisenberg_model", "nn_cutoff"],
            )
        )
        hmodels = [hmodel for hmodel in docs["heisenberg_model"]]
        cutoffs = [cutoff for cutoff in docs["nn_cutoff"]]

        # Check for J_ij convergence
        converged_list = []
        for cutoff, hmodel in zip(cutoffs, hmodels):
            ex_params = hmodel.ex_params
            converged = True

            # J_ij exists
            if len(ex_params) > 1:
                E0 = abs(ex_params["E0"])
                for k, v in ex_params.items():
                    if k != "E0":
                        if abs(v) > E0:  # |J_ij| > |E0| unphysical!
                            converged = False
            else:  # Only <J> was computed
                converged = False

            if converged:
                converged_list.append(hmodel)

        # If multiple Heisenberg Models converged, take the maximal cutoff
        if len(converged_list) > 0:
            hmodel = converged_list[-1]  # Largest cutoff
        else:
            hmodel = None

        # Update FW spec with converged hmodel or None
        update_spec = {"converged_heisenberg_model": hmodel}

        return FWAction(update_spec=update_spec)


@explicit_serialize
class VampireMC(FiretaskBase):
    """
    If there is a suitable Heisenberg model, run Monte Carlo to compute the critical temperature.

    Args:
        db_file (str): path to file containing the database credentials.
        exchange_wf_uuid (int): Unique id for record keeping.

    """

    required_params = ["db_file", "exchange_wf_uuid", "parent_structure"]
    optional_params = []

    def run_task(self, fw_spec):

        db_file = env_chk(self["db_file"], fw_spec)
        wf_uuid = self["exchange_wf_uuid"]
        formula_pretty = self["parent_structure"].composition.reduced_formula

        # Get Heisenberg models from db
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        mmdb.collection = mmdb.db["exchange"]

        task_doc = {"wf_meta.wf_uuid": wf_uuid, "formula_pretty": formula_pretty}

        # Get a converged Heisenberg model if one was found
        hmodel = fw_spec["converged_heisenberg_model"]

        if hmodel:
            vc = VampireCaller(hm=hmodel)
            vo = vc.output
            task_doc["vampire_output"] = vo
        else:
            task_doc["vampire_output"] = None

        # Insert vampire output into exchange collection
        mmdb.collection.insert(task_doc)
