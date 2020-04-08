from fireworks import FiretaskBase, FWAction, explicit_serialize

from monty.serialization import loadfn, dumpfn

from atomate.utils.utils import get_logger, env_chk
from atomate.vasp.database import VaspCalcDb

from datetime import datetime
import numpy as np

from pymatgen import Structure
from pymatgen.command_line.bader_caller import bader_analysis_from_path
from pymatgen.analysis.magnetism.heisenberg import HeisenbergMapper, HeisenbergModel
from pymatgen.analysis.magnetism import CollinearMagneticStructureAnalyzer, Ordering
from pymatgen.command_line.vampire_caller import VampireCaller

__author__ = "Nathan C. Frey"
__email__ = "ncfrey@lbl.gov"

logger = get_logger(__name__)


@explicit_serialize
class HeisenbergModelMapping(FiretaskBase):
    """
    Map structures and energies to a Heisenberg model and compute exchange
    parameters for a given NN cutoff.

    * heisenberg_settings: 
        cutoff (float): Starting point for nearest neighbor search.
        tol (float): Tolerance for equivalent NN bonds.

    Args:
        structures (list): Magnetic structures.
        energies (list): Energies / atom (eV).

    Optional parameters:
        heisenberg_settings (dict): A config dict for Heisenberg model 
            mapping, detailed above.

    """

    required_params = ["structures", "energies"]

    optional_params = ["heisenberg_settings"]

    def run_task(self, fw_spec):

        structures = self["structures"]
        energies = self["energies"]

        heisenberg_settings = self.get("heisenberg_settings", {})

        # Total energies
        energies = [e * len(s) for e, s in zip(energies, structures)]

        # Map system to a Heisenberg Model
        hmapper = HeisenbergMapper(structures, energies, **heisenberg_settings)

        # Get MSONable Heisenberg Model
        hmodel = hmapper.get_heisenberg_model()
        hmodel_dict = hmodel.as_dict()

        # Update FW spec with model
        name = "heisenberg_model_" + str(hmodel.cutoff).replace(".", "_")
        update_spec = {name: hmodel}

        # Write to file
        dumpfn(hmodel_dict, "heisenberg_model.json")

        return FWAction(update_spec=update_spec)


@explicit_serialize
class HeisenbergModelToDb(FiretaskBase):
    """
    Insert Heisenberg Model object into a DB. Assumes you are in a
    directory with a model written to a .json file.

    Args:
        db_file (str): path to file containing the database credentials.
        wf_uuid (int): Unique id for record keeping.

    """

    required_params = ["db_file", "wf_uuid"]

    def run_task(self, fw_spec):

        db_file = env_chk(self["db_file"], fw_spec)
        wf_uuid = self["wf_uuid"]

        hmodel = loadfn("heisenberg_model.json")

        # Exchange collection
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        mmdb.collection = mmdb.db["exchange"]

        hmodel_dict = hmodel.as_dict()

        parent_structure = hmodel.structures[0]
        formula_pretty = parent_structure.composition.reduced_formula

        wf_meta = {"wf_uuid": wf_uuid}
        task_doc = {
            "wf_meta": wf_meta,
            "formula_pretty": formula_pretty,
            "nn_cutoff": hmodel.cutoff,
            "nn_tol": hmodel.tol,
            "heisenberg_model": hmodel_dict,
            "task_name": "heisenberg model",
        }

        if fw_spec.get("tags", None):
            task_doc["tags"] = fw_spec["tags"]

        mmdb.collection.insert_one(task_doc)


@explicit_serialize
class HeisenbergConvergence(FiretaskBase):
    """
    Quick check to see if the Heisenberg model has "converged" for any
    particular nearest neighbor cutoff value in the sense that |J_ij| < |E_0|
    for all i, j. 

    If not, it doesn't make sense to do Monte Carlo and the system is either
    1) not well described by the Heisenberg Model (perhaps the magnetic
    moments are not localized) or 2) not correctly captured by the provided
    low-energy magnetic orderings.

    Args:
        db_file (str): path to file containing the database credentials.
        wf_uuid (int): Unique id for record keeping.
        average (bool): <J> exchange param only.

    TODO:
        * More robust convergence check

    """

    required_params = ["db_file", "wf_uuid", "average"]

    def run_task(self, fw_spec):

        db_file = env_chk(self["db_file"], fw_spec)
        wf_uuid = self["wf_uuid"]

        # Get Heisenberg models from db
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        mmdb.collection = mmdb.db["exchange"]

        # Get documents
        docs = list(
            mmdb.collection.find(
                {"wf_meta.wf_uuid": wf_uuid}, ["heisenberg_model", "nn_cutoff"]
            )
        )

        hmodels = [HeisenbergModel.from_dict(d["heisenberg_model"]) for d in docs]
        cutoffs = [d["nn_cutoff"] for d in docs]

        # Take <J> only if no static calcs
        if self["average"]:
            hmodel = hmodels[0].as_dict()
        else:
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
                    converged_list.append(hmodel.as_dict())

            # If multiple Heisenberg Models converged, take the maximal cutoff
            if len(converged_list) > 0:
                hmodel = converged_list[-1]  # Largest cutoff
            else:  # No converged models
                hmodel = None

        # Update FW spec with converged hmodel or None
        update_spec = {"converged_heisenberg_model": hmodel}

        return FWAction(update_spec=update_spec)


@explicit_serialize
class VampireMC(FiretaskBase):
    """
    If there is a suitable Heisenberg model, run Monte Carlo to compute the
    critical temperature.

    Args:
        db_file (str): path to file containing the database credentials.
        wf_uuid (int): Unique id for record keeping.

    Optional args:
        mc_settings (dict): A configuration dict for monte carlo.
            See pymatgen.command_line.VampireCaller for options.

    TODO:
        * Include HeisenbergModel convergence check.

    """

    required_params = ["db_file", "wf_uuid"]

    optional_params = ["mc_settings"]

    def run_task(self, fw_spec):

        db_file = env_chk(self["db_file"], fw_spec)
        wf_uuid = self["wf_uuid"]
        mc_settings = self.get("mc_settings", {})

        # Get Heisenberg models from db
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        mmdb.collection = mmdb.db["exchange"]

        # Get documents
        docs = list(
            mmdb.collection.find(
                {"wf_meta.wf_uuid": wf_uuid}, ["heisenberg_model", "nn_cutoff"]
            )
        )

        hmodels = [HeisenbergModel.from_dict(d["heisenberg_model"]) for d in docs]
        cutoffs = [hmodel.cutoff for hmodel in hmodels]
        ordered_hmodels = [h for _, h in sorted(zip(cutoffs, hmodels), reverse=False)]
        # Take the model with smallest NN cutoff
        hmodel = ordered_hmodels[0]

        # Get a converged Heisenberg model if one was found
        # if fw_spec["converged_heisenberg_model"]:
        #     hmodel = HeisenbergModel.from_dict(fw_spec["converged_heisenberg_model"])

        vc = VampireCaller(hm=hmodel, **mc_settings)
        vampire_output = vc.output

        # Update FW spec
        update_spec = {"vampire_output": vampire_output}

        # Write to file
        dumpfn(vampire_output.as_dict(), "vampire_output.json")

        return FWAction(update_spec=update_spec)


@explicit_serialize
class VampireToDb(FiretaskBase):
    """Insert VAMPIRE output into DB.

    Assumes you are in a directory with output written to a .json.

    Args:
        db_file (str): path to file containing the database credentials.
        wf_uuid (int): Unique id for record keeping.

    """

    required_params = ["db_file", "wf_uuid"]

    def run_task(self, fw_spec):

        db_file = env_chk(self["db_file"], fw_spec)
        wf_uuid = self["wf_uuid"]
        vampire_output = loadfn("vampire_output.json")

        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        mmdb.collection = mmdb.db["exchange"]

        task_doc = {"wf_meta": {"wf_uuid": wf_uuid}, "task_name": "vampire caller"}

        if fw_spec.get("tags", None):
            task_doc["tags"] = fw_spec["tags"]

        task_doc["vampire_output"] = vampire_output.as_dict()
        mmdb.collection.insert_one(task_doc)
