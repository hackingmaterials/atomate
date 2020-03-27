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

    Args:
        db_file (str): path to file containing the database credentials.
        wf_uuid (int): Unique id for record keeping.
        cutoff (float): NN search cutoff (Angstrom).
        tol (float): distance tolerance for similar bonds.
        average (bool): <J> only or full NN, NNN HeisenbergModel.

    Optional parameters:
        structures (list): Magnetic structures.
        energies (list): Energies / atom (eV).

    """

    required_params = [
        "db_file",
        "wf_uuid",
        "cutoff",
        "tol",
        "average",
    ]

    # If only doing <J>, give the original structure/energy inputs
    optional_params = ["structures", "energies"]

    def run_task(self, fw_spec):

        db_file = env_chk(self["db_file"], fw_spec)
        wf_uuid = self["wf_uuid"]

        # Get magnetic orderings collection from db
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        mmdb.collection = mmdb.db["magnetic_orderings"]

        # Look up static calcs if needed
        if not self["average"]:

            # Get documents
            docs = list(
                mmdb.collection.find(
                    {"wf_meta.wf_uuid": wf_uuid},
                    ["task_id", "structure", "energy_per_atom"],
                )
            )

            # Get structures and energy / unit cell
            structures = [Structure.from_dict(d["structure"]) for d in docs]
            epas = [d["energy_per_atom"] for d in docs]
            n_atoms = len(structures[0])
            energies = [e * n_atoms for e in epas]
        else:
            structures = self["structures"]
            energies = self["energies"]

        # Total energies
        energies = [e * len(s) for e, s in zip(energies, structures)]

        # Map system to a Heisenberg Model
        hmapper = HeisenbergMapper(structures, energies, self["cutoff"], self["tol"])

        # Get MSONable Heisenberg Model
        hmodel = hmapper.get_heisenberg_model()
        hmodel_dict = hmodel.as_dict()

        # Update FW spec with model
        name = "heisenberg_model_" + str(self["cutoff"])
        update_spec = {name: hmodel}

        # Write to file
        dumpfn(hmodel_dict, name + ".json")

        return FWAction(update_spec=update_spec)


@explicit_serialize
class HeisenbergModelToDb(FiretaskBase):
    """
    Insert Heisenberg Model object into a DB. Assumes you are in a
    directory with a model written to a .json file.

    Args:
        db_file (str): path to file containing the database credentials.
        wf_uuid (int): Unique id for record keeping.
        cutoff (float): NN search cutoff (Angstrom).

    """

    required_params = [
        "db_file",
        "wf_uuid",
        "cutoff"
    ]

    def run_task(self, fw_spec):

        db_file = env_chk(self["db_file"], fw_spec)
        wf_uuid = self["wf_uuid"]

        name = "heisenberg_model_" + str(self["cutoff"])
        hmodel = loadfn(name + ".json")
        #hmodel = HeisenbergModel.from_dict(hmodel_dict)

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
            "nn_cutoff": self["cutoff"],
            "nn_tol": hmodel.tol,
            "heisenberg_model": hmodel_dict,
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
        mc_settings (dict): A configuration dict for monte carlo.
        average (bool): Only <J> exchange param.

    TODO:
        * Include HeisenbergModel convergence check.

    """

    required_params = [
        "db_file",
        "wf_uuid",
        "mc_settings",
        "average",
    ]

    def run_task(self, fw_spec):

        db_file = env_chk(self["db_file"], fw_spec)
        wf_uuid = self["wf_uuid"]
        mc_settings = self["mc_settings"]
        average = self["average"]

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
        hmodel = hmodels[0]  # Take the model with smallest NN cutoff

        # Vampire monte carlo settings
        mc_box_size = mc_settings["mc_box_size"]
        equil_timesteps = mc_settings["equil_timesteps"]
        mc_timesteps = mc_settings["mc_timesteps"]

        # Get a converged Heisenberg model if one was found
        # if fw_spec["converged_heisenberg_model"]:
        #     hmodel = HeisenbergModel.from_dict(fw_spec["converged_heisenberg_model"])

        vc = VampireCaller(
            mc_box_size=mc_box_size,
            equil_timesteps=equil_timesteps,
            mc_timesteps=mc_timesteps,
            hm=hmodel,
            avg=average,
        )
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

    required_params = [
        "db_file",
        "wf_uuid",
    ]

    def run_task(self, fw_spec):

        db_file = env_chk(self["db_file"], fw_spec)
        wf_uuid = self["wf_uuid"]
        vampire_output = loadfn("vampire_output.json")

        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        mmdb.collection = mmdb.db["exchange"]

        task_doc = {"wf_meta": {"wf_uuid": wf_uuid}}

        if fw_spec.get("tags", None):
            task_doc["tags"] = fw_spec["tags"]

        task_doc["vampire_output"] = vampire_output.as_dict()
        mmdb.collection.insert_one(task_doc)