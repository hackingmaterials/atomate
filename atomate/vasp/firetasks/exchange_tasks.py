from fireworks import FiretaskBase, FWAction, explicit_serialize

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
class ExchangeToDB(FiretaskBase):
    """
    Used to aggregate tasks docs from magnetic ordering workflow.

    Args:
        db_file (str): path to the db file that holds your tasks
        collection and that you want to hold the magnetic_orderings
        collection
        wf_uuid (str): auto-generated from get_wf_magnetic_orderings,
        used to make it easier to retrieve task docs
        parent_structure: Structure of parent crystal (not magnetically
        ordered)
    """

    required_params = [
        "db_file",
        "wf_uuid",
        "parent_structure",
        "perform_bader",
        "scan",
    ]
    optional_params = ["origins", "input_index"]

    def run_task(self, fw_spec):

        uuid = self["wf_uuid"]
        db_file = env_chk(self.get("db_file"), fw_spec)
        to_db = self.get("to_db", True)

        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)

        formula = self["parent_structure"].formula
        formula_pretty = self["parent_structure"].composition.reduced_formula

        # get ground state energy
        task_label_regex = "static" if not self["scan"] else "optimize"
        docs = list(
            mmdb.collection.find(
                {"wf_meta.wf_uuid": uuid, "task_label": {"$regex": task_label_regex}},
                ["task_id", "output.energy_per_atom"],
            )
        )

        energies = [d["output"]["energy_per_atom"] for d in docs]
        ground_state_energy = min(energies)
        idx = energies.index(ground_state_energy)
        ground_state_task_id = docs[idx]["task_id"]
        if energies.count(ground_state_energy) > 1:
            logger.warning(
                "Multiple identical energies exist, "
                "duplicate calculations for {}?".format(formula)
            )

        # get results for different orderings
        docs = list(
            mmdb.collection.find(
                {"task_label": {"$regex": task_label_regex}, "wf_meta.wf_uuid": uuid}
            )
        )

        summaries = []

        for d in docs:

            optimize_task_label = d["task_label"]
            optimize_task = dict(
                mmdb.collection.find_one(
                    {"wf_meta.wf_uuid": uuid, "task_label": optimize_task_label}
                )
            )
            input_structure = Structure.from_dict(optimize_task["input"]["structure"])
            input_magmoms = optimize_task["input"]["incar"]["MAGMOM"]
            input_structure.add_site_property("magmom", input_magmoms)

            final_structure = Structure.from_dict(d["output"]["structure"])

            input_analyzer = CollinearMagneticStructureAnalyzer(
                input_structure, threshold=0.0
            )
            final_analyzer = CollinearMagneticStructureAnalyzer(
                final_structure, threshold=0.0
            )

            if d["task_id"] == ground_state_task_id:
                stable = True
                decomposes_to = None
            else:
                stable = False
                decomposes_to = ground_state_task_id
            energy_above_ground_state_per_atom = (
                d["output"]["energy_per_atom"] - ground_state_energy
            )
            energy_diff_relax_static = (
                optimize_task["output"]["energy_per_atom"]
                - d["output"]["energy_per_atom"]
            )

            # tells us the order in which structure was guessed
            # 1 is FM, then AFM..., -1 means it was entered manually
            # useful to give us statistics about how many orderings
            # we actually need to calculate
            task_label = d["task_label"].split(" ")
            ordering_index = task_label.index("ordering")
            ordering_index = int(task_label[ordering_index + 1])
            if self.get("origins", None):
                ordering_origin = self["origins"][ordering_index]
            else:
                ordering_origin = None

            final_magmoms = final_structure.site_properties["magmom"]
            magmoms = {"vasp": final_magmoms}
            if self["perform_bader"]:
                # if bader has already been run during task ingestion,
                # use existing analysis
                if "bader" in d:
                    magmoms["bader"] = d["bader"]["magmom"]
                # else try to run it
                else:
                    try:
                        dir_name = d["dir_name"]
                        # strip hostname if present, implicitly assumes
                        # ToDB task has access to appropriate dir
                        if ":" in dir_name:
                            dir_name = dir_name.split(":")[1]
                        magmoms["bader"] = bader_analysis_from_path(dir_name)["magmom"]
                        # prefer bader magmoms if we have them
                        final_magmoms = magmoms["bader"]
                    except Exception as e:
                        magmoms["bader"] = "Bader analysis failed: {}".format(e)

            # Specify very small threshold
            input_order_check = [0 if abs(m) < 0.001 else m for m in input_magmoms]
            final_order_check = [0 if abs(m) < 0.001 else m for m in final_magmoms]
            ordering_changed = not np.array_equal(
                np.sign(input_order_check), np.sign(final_order_check)
            )

            symmetry_changed = (
                final_structure.get_space_group_info()[0]
                != input_structure.get_space_group_info()[0]
            )

            total_magnetization = abs(
                d["calcs_reversed"][0]["output"]["outcar"]["total_magnetization"]
            )
            num_formula_units = sum(
                d["calcs_reversed"][0]["composition_reduced"].values()
            ) / sum(d["calcs_reversed"][0]["composition_unit_cell"].values())
            total_magnetization_per_formula_unit = (
                total_magnetization / num_formula_units
            )
            total_magnetization_per_unit_volume = (
                total_magnetization / final_structure.volume
            )

            summary = {
                "formula": formula,
                "formula_pretty": formula_pretty,
                "parent_structure": self["parent_structure"].as_dict(),
                "wf_meta": d["wf_meta"],  # book-keeping
                "task_id": d["task_id"],
                "structure": final_structure.as_dict(),
                "magmoms": magmoms,
                "input": {
                    "structure": input_structure.as_dict(),
                    "ordering": input_analyzer.ordering.value,
                    "symmetry": input_structure.get_space_group_info()[0],
                    "index": ordering_index,
                    "origin": ordering_origin,
                    "input_index": self.get("input_index", None),
                },
                "total_magnetization": total_magnetization,
                "total_magnetization_per_formula_unit": total_magnetization_per_formula_unit,
                "total_magnetization_per_unit_volume": total_magnetization_per_unit_volume,
                "ordering": final_analyzer.ordering.value,
                "ordering_changed": ordering_changed,
                "symmetry": final_structure.get_space_group_info()[0],
                "symmetry_changed": symmetry_changed,
                "energy_per_atom": d["output"]["energy_per_atom"],
                "stable": stable,
                "decomposes_to": decomposes_to,
                "energy_above_ground_state_per_atom": energy_above_ground_state_per_atom,
                "energy_diff_relax_static": energy_diff_relax_static,
                "created_at": datetime.utcnow(),
            }

            if fw_spec.get("tags", None):
                summary["tags"] = fw_spec["tags"]

            summaries.append(summary)

        mmdb.collection = mmdb.db["magnetic_orderings"]
        mmdb.collection.insert(summaries)

        logger.info("Magnetic orderings calculation complete.")


@explicit_serialize
class HeisenbergModelMapping(FiretaskBase):
    """
    Map structures and energies to a Heisenberg model and compute exchange parameters for a given NN cutoff.

    Args:
        db_file (str): path to file containing the database credentials.
        exchange_wf_uuid (int): Unique id for record keeping.
        cutoff (float): NN search cutoff (Angstrom).
        tol (float): distance tolerance for similar bonds.
        avg (bool): <J> only or full NN, NNN HeisenbergModel.
        structures (list): Magnetic structures.
        energies (list): Energies / atom (eV).
    """

    required_params = [
        "db_file",
        "exchange_wf_uuid",
        "parent_structure",
        "cutoff",
        "tol",
        "avg",
    ]

    # If only doing <J>, give the original structure/energy inputs
    optional_params = [
        "structures",
        "energies",
        ]

    def run_task(self, fw_spec):

        db_file = env_chk(self["db_file"], fw_spec)
        wf_uuid = self["exchange_wf_uuid"]
        formula = self["parent_structure"].formula
        formula_pretty = self["parent_structure"].composition.reduced_formula

        # Get magnetic orderings collection from db
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        mmdb.collection = mmdb.db["magnetic_orderings"]

        # Look up static calcs if needed
        if not self["avg"]:

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
        name = "heisenberg_model_" + str(self["cutoff"])

        wf_meta = {"wf_uuid": wf_uuid}
        task_doc = {
            "wf_meta": wf_meta,
            "formula_pretty": formula_pretty,
            "nn_cutoff": self["cutoff"],
            "nn_tol": self["tol"],
            "heisenberg_model": hmodel_dict,
        }

        if fw_spec.get("tags", None):
            task_doc["tags"] = fw_spec["tags"]

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
        parent_structure (Structure): Structure object.
        avg (bool): <J> exchange param only.

    TODO:
        * More robust convergence check

    """

    required_params = ["db_file", "exchange_wf_uuid", "parent_structure", "avg"]
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
                {"wf_meta.wf_uuid": wf_uuid}, ["heisenberg_model", "nn_cutoff"]
            )
        )

        hmodels = [HeisenbergModel.from_dict(d["heisenberg_model"]) for d in docs]
        cutoffs = [d["nn_cutoff"] for d in docs]

        # Take <J> only if no static calcs
        if self["avg"]:
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
    If there is a suitable Heisenberg model, run Monte Carlo to compute the critical temperature.

    Args:
        db_file (str): path to file containing the database credentials.
        exchange_wf_uuid (int): Unique id for record keeping.
        parent_structure (Structure): Structure object with magmoms.
        mc_settings (dict): A configuration dict for monte carlo.
        avg (bool): Only <J> exchange param.

    TODO:
        * Include HeisenbergModel convergence check.

    """

    required_params = ["db_file", "exchange_wf_uuid", "parent_structure", "mc_settings", "avg"]
    optional_params = []

    def run_task(self, fw_spec):

        db_file = env_chk(self["db_file"], fw_spec)
        wf_uuid = self["exchange_wf_uuid"]
        formula_pretty = self["parent_structure"].composition.reduced_formula
        mc_settings = self["mc_settings"]
        avg = self["avg"]

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

        task_doc = {"wf_meta": {"wf_uuid": wf_uuid}, "formula_pretty": formula_pretty}

        if fw_spec.get("tags", None):
            task_doc["tags"] = fw_spec["tags"]

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
            avg=avg,
        )
        vo = vc.output
        task_doc["vampire_output"] = vo.as_dict()
        
        # else:
        #     task_doc["vampire_output"] = None

        # Insert vampire output into exchange collection
        mmdb.collection.insert(task_doc)
