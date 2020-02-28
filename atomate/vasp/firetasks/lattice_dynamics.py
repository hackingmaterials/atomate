import json
import os
import subprocess
import shlex
from datetime import datetime
from pathlib import Path

import numpy as np
from monty.dev import requires
from monty.serialization import dumpfn, loadfn
from pymongo import ReturnDocument

from atomate.utils.utils import env_chk, get_logger
from atomate.vasp.analysis.lattice_dynamics import (
    FIT_METHOD,
    IMAGINARY_TOL,
    MAX_IMAGINARY_FREQ,
    MAX_N_IMAGINARY,
    fit_force_constants,
    get_cutoffs,
)
from atomate.vasp.database import VaspCalcDb
from fireworks import FiretaskBase, FWAction, explicit_serialize
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.phonopy import get_phonon_band_structure_from_fc, \
    get_phonon_dos_from_fc, get_phonon_band_structure_symm_line_from_fc
from pymatgen.io.shengbte import Control
from pymatgen.transformations.standard_transformations import (
    SupercellTransformation,
)

try:
    import hiphive
except ImportError:
    hiphive = False


__author__ = "Alex Ganose"
__email__ = "aganose@lbl.gov"

logger = get_logger(__name__)

# define shared constants
DEFAULT_TEMPERATURE = 300
MESH_DENSITY = 100.0  # should always be a float


@explicit_serialize
class CollectPerturbedStructures(FiretaskBase):
    """
    Aggregate the structures and forces of perturbed supercells.

    Requires a ``perturbed_tasks`` key to be set in the fw_spec, containing a
    list of dictionaries, each with the keys:

    - "parent_structure": The original structure.
    - "supercell_matrix": The supercell transformation matrix.
    - "structure": The perturbed supercell structure.
    - "forces": The forces on each site in the supercell.
    """

    def run_task(self, fw_spec):
        results = fw_spec.get("perturbed_tasks", [])

        if len(results) == 0:
            # can happen if all parents fizzled
            raise RuntimeError("No perturbed tasks found in firework spec")

        logger.info("Found {} perturbed structures".format(len(results)))

        structure = results[0]["parent_structure"]
        supercell_matrix = results[0]["supercell_matrix"]

        structures = [r["structure"] for r in results]
        forces = [np.asarray(r["forces"]) for r in results]

        # regenerate pure supercell structure
        st = SupercellTransformation(supercell_matrix)
        supercell_structure = st.apply_transformation(structure)

        structure_data = {
            "structure": structure,
            "supercell_structure": supercell_structure,
            "supercell_matrix": supercell_matrix,
        }

        logger.info("Writing structure and forces data.")
        dumpfn(structures, "perturbed_structures.json")
        dumpfn(forces, "perturbed_forces.json")
        dumpfn(structure_data, "structure_data.json")


@explicit_serialize
class RunHiPhive(FiretaskBase):
    """
    Fit force constants using hiPhive.

    Requires "perturbed_structures.json", "perturbed_forces.json", and
    "structure_data.json" files to be present in the current working directory.

    Optional parameters:
        cutoffs (Optional[list[list]]): A list of cutoffs to trial. If None,
            a set of trial cutoffs will be generated based on the structure
            (default).
        imaginary_tol (float): Tolerance used to decide if a phonon mode
            is imaginary, in THz.
        max_n_imaginary (int): Maximum number of imaginary modes allowed in the
            the final fitted force constant solution. If this criteria is not
            reached by any cutoff combination this FireTask will fizzle.
        max_imaginary_freq (float): Maximum allowed imaginary frequency in the
            final fitted force constant solution. If this criteria is not
            reached by any cutoff combination this FireTask will fizzle.
        fit_method (str): Method used for fitting force constants. This can be
            any of the values allowed by the hiphive ``Optimizer`` class.
    """

    optional_params = [
        "cutoffs",
        "imaginary_tol",
        "max_n_imaginary",
        "max_imaginary_freq",
        "fit_method",
    ]

    @requires(hiphive, "hiphive is required for lattice dynamics workflow")
    def run_task(self, fw_spec):
        from hiphive.utilities import get_displacements

        max_n_imaginary = self.get("max_n_imaginary", MAX_N_IMAGINARY)
        max_imaginary_freq = self.get("max_imaginary_freq", MAX_IMAGINARY_FREQ)
        imaginary_tol = self.get("imaginary_tol", IMAGINARY_TOL)
        fit_method = self.get("fit_method", FIT_METHOD)

        all_structures = loadfn("perturbed_structures.json")
        all_forces = loadfn("perturbed_forces.json")
        structure_data = loadfn("structure_data.json")
        parent_structure = structure_data["structure"]
        supercell_matrix = structure_data["supercell_matrix"]
        supercell_structure = structure_data["supercell_structure"]

        structures = []
        supercell_atoms = AseAtomsAdaptor.get_atoms(supercell_structure)
        for structure, forces in zip(all_structures, all_forces):
            atoms = AseAtomsAdaptor.get_atoms(structure)
            displacements = get_displacements(atoms, supercell_atoms)
            atoms.new_array("displacements", displacements)
            atoms.new_array("forces", forces)
            atoms.positions = supercell_atoms.get_positions()
            structures.append(atoms)

        cutoffs = self.get("cutoffs") or get_cutoffs(supercell_structure)
        force_constants, fitting_data = fit_force_constants(
            parent_structure,
            supercell_matrix,
            structures,
            cutoffs,
            imaginary_tol,
            max_n_imaginary,
            max_imaginary_freq,
            fit_method,
        )

        dumpfn(fitting_data, "fitting_data.json")

        if force_constants is None:
            # fitting failed
            raise RuntimeError(
                "Could not find a force constant solution with less than {} "
                "imaginary modes.\n"
                "Fitting results: {}".format(max_n_imaginary, fitting_data)
            )

        else:
            logger.info("Writing force constants.")
            force_constants.write("force_constants.fcs")

            atoms = AseAtomsAdaptor.get_atoms(parent_structure)
            force_constants.write_to_shengBTE("FORCE_CONSTANTS_3RD", atoms)
            force_constants.write_to_phonopy(
                "FORCE_CONSTANTS_2ND", format="text"
            )


@explicit_serialize
class ForceConstantsToDb(FiretaskBase):
    """
    Add force constants, phonon band structure and density of states
    to the database.

    Assumes you are in a directory with the force constants, fitting
    data, and structure data written to files.

    Required parameters:
        db_file (str): Path to DB file for the database that contains the
            perturbed structure calculations.

    Optional parameters:
        mesh_density (float): The density of the q-point mesh used to calculate
            the phonon density of states. See the docstring for the ``mesh``
            argument in Phonopy.init_mesh() for more details.
        additional_fields (dict): Additional fields added to the document, such
            as user-defined tags, name, ids, etc.
    """

    required_params = ["db_file"]
    optional_params = ["mesh_density", "additional_fields"]

    @requires(hiphive, "hiphive is required for lattice dynamics workflow")
    def run_task(self, fw_spec):
        from hiphive.force_constants import ForceConstants

        db_file = env_chk(self.get("db_file"), fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)

        fc = ForceConstants.read("force_constants.fcs")
        fitting_data = loadfn("fitting_data.json")
        structure_data = loadfn("structure_data.json")
        forces = loadfn("perturbed_forces.json")
        structures = loadfn("perturbed_structures.json")

        structure = structure_data["structure"]
        supercell_structure = structure_data["supercell_structure"]
        supercell_matrix = structure_data["supercell_matrix"]

        phonopy_fc = fc.get_fc_array(order=2)

        logger.info("Getting uniform phonon band structure.")
        uniform_bs = get_phonon_band_structure_from_fc(
            structure, supercell_matrix, phonopy_fc
        )

        logger.info("Getting line mode phonon band structure.")
        lm_bs = get_phonon_band_structure_symm_line_from_fc(
            structure, supercell_matrix, phonopy_fc
        )

        logger.info("Getting phonon density of states.")
        mesh_density = self.get("mesh_density", MESH_DENSITY)
        dos = get_phonon_dos_from_fc(
            structure, supercell_matrix, phonopy_fc, mesh_density=mesh_density
        )

        logger.info("Inserting phonon objects into database.")
        dos_fsid, _ = mmdb.insert_gridfs(
            dos.to_json(), collection="phonon_dos_fs"
        )
        uniform_bs_fsid, _ = mmdb.insert_gridfs(
            uniform_bs.to_json(), collection="phonon_bandstructure_fs"
        )
        lm_bs_fsid, _ = mmdb.insert_gridfs(
            lm_bs.to_json(), collection="phonon_bandstructure_fs"
        )

        logger.info("Inserting force constants into database.")
        fc_json = json.dumps(
            {str(k): v.tolist() for k, v in fc.get_fc_dict().items()}
        )
        fc_fsid, _ = mmdb.insert_gridfs(
            fc_json, collection="phonon_force_constants_fs"
        )

        data = {
            "structure": structure.as_dict(),
            "supercell_matrix": supercell_matrix,
            "supercell_structure": supercell_structure.as_dict(),
            "perturbed_structures": [s.as_dict() for s in structures],
            "perturbed_forces": [f.tolist() for f in forces],
            "fitting_data": fitting_data,
            "force_constants_fs_id": fc_fsid,
            "tags": fw_spec.get("tags", None),
            "formula_pretty": structure.composition.reduced_formula,
            "phonon_dos_fs_id": dos_fsid,
            "phonon_bandstructure_uniform_fs_id": uniform_bs_fsid,
            "phonon_bandstructure_line_fs_id": lm_bs_fsid,
            "created_at": datetime.utcnow(),
        }
        data.update(self.get("additional_fields", {}))

        # Get an id for the force constants
        fitting_id = _get_fc_fitting_id(mmdb)
        metadata = {"fc_fitting_id": fitting_id, "fc_fitting_dir": os.getcwd()}
        data.update(metadata)

        mmdb.db.lattice_dynamics.insert_one(data)
        logger.info("Finished inserting force constants")

        return FWAction(update_spec=metadata)


@explicit_serialize
class RunShengBTE(FiretaskBase):
    """
    Run ShengBTE to calculate lattice thermal conductivity. Presumes
    the FORCE_CONSTANTS_3RD and FORCE_CONSTANTS_2ND, and a "structure_data.json"
    file, with the keys "structure", " and "supercell_matrix" is in the current
    directory.

    Required parameters:
        shengbte_cmd (str): The name of the shengbte executable to run. Supports
            env_chk.

    Optional parameters:
        temperature (float or dict): The temperature to calculate the lattice
            thermal conductivity for. Can be given as a single float, or a
            dictionary with the keys "min", "max", "step".
        control_kwargs (dict): Options to be included in the ShengBTE control
            file.
    """

    required_params = ["shengbte_cmd"]
    optional_params = ["temperature", "control_kwargs"]

    def run_task(self, fw_spec):
        structure_data = loadfn("structure_data.json")
        structure = structure_data["structure"]
        supercell_matrix = structure_data["supercell_matrix"]

        control_dict = {
            "scalebroad": 0.5,
            "nonanalytic": False,
            "isotopes": False,
            "temperature": self.get("temperature", DEFAULT_TEMPERATURE),
            "scell": np.diag(supercell_matrix).tolist(),
        }
        control_kwargs = self.get("control_kwargs") or {}
        control_dict.update(control_kwargs)
        control = Control().from_structure(structure, **control_dict)
        control.to_file()

        shengbte_cmd = env_chk(self["shengbte_cmd"], fw_spec)

        if isinstance(shengbte_cmd, str):
            shengbte_cmd = os.path.expandvars(shengbte_cmd)
            shengbte_cmd = shlex.split(shengbte_cmd)

        shengbte_cmd = list(shengbte_cmd)
        logger.info("Running command: {}".format(shengbte_cmd))

        with open("shengbte.out", "w") as f_std, open(
            "shengbte_err.txt", "w", buffering=1
        ) as f_err:
            # use line buffering for stderr
            return_code = subprocess.call(
                shengbte_cmd, stdout=f_std, stderr=f_err
            )
        logger.info(
            "Command {} finished running with returncode: {}".format(
                shengbte_cmd, return_code
            )
        )

        if return_code == 1:
            raise RuntimeError(
                "Running ShengBTE failed. Check '{}/shengbte_err.txt' for "
                "details.".format(os.getcwd())
            )


@explicit_serialize
class ShengBTEToDb(FiretaskBase):
    """
    Add lattice thermal conductivity results to database.

    Assumes you are in a directory with the ShengBTE results in.

    Required parameters:
        db_file (str): Path to DB file for the database that contains the
            perturbed structure calculations.

    Optional parameters:
        additional_fields (dict): Additional fields added to the document.
    """

    required_params = ["db_file"]
    optional_params = ["additional_fields"]

    def run_task(self, fw_spec):
        db_file = env_chk(self.get("db_file"), fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)

        control = Control().from_file("CONTROL")
        structure = control.get_structure()
        supercell_matrix = np.diag(control["scell"])

        if Path("BTE.KappaTensorVsT_CONV").exists():
            filename = "BTE.KappaTensorVsT_CONV"
        elif Path("BTE.KappaTensorVsT_RTA").exists():
            filename = "BTE.KappaTensorVsT_RTA"
        else:
            raise RuntimeError("Could not find ShengBTE output files.")

        bte_data = np.loadtxt(filename)
        if len(bte_data.shape) == 1:
            # pad extra axis to make compatible with multiple temperatures
            bte_data = bte_data[None, :]

        temperatures = bte_data[:, 0].tolist()
        kappa = bte_data[:, 1:10].reshape(-1, 3, 3).tolist()

        data = {
            "structure": structure.as_dict(),
            "supercell_matrix": supercell_matrix.tolist(),
            "temperatures": temperatures,
            "lattice_thermal_conductivity": kappa,
            "control": control.as_dict(),
            "tags": fw_spec.get("tags", None),
            "formula_pretty": structure.composition.reduced_formula,
            "shengbte_dir": os.getcwd(),
            "fc_fitting_id": fw_spec.get("fc_fitting_id", None),
            "fc_fitting_dir": fw_spec.get("fc_fitting_dir", None),
            "created_at": datetime.utcnow(),
        }
        data.update(self.get("additional_fields", {}))

        mmdb.collection = mmdb.db["lattice_thermal_conductivity"]
        mmdb.collection.insert(data)


def _get_fc_fitting_id(mmdb: VaspCalcDb) -> int:
    """Helper method to get a force constant fitting id."""
    fc_id = mmdb.db.counter.find_one_and_update(
        {"_id": "fc_fitting_id"},
        {"$inc": {"c": 1}},
        return_document=ReturnDocument.AFTER,
    )
    if fc_id is None:
        mmdb.db.counter.insert({"_id": "fc_fitting_id", "c": 1})
        fc_id = 1
    else:
        fc_id = fc_id["c"]

    return fc_id
