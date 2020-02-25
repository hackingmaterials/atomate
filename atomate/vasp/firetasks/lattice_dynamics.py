import os
import shlex

from pathlib import Path

from uuid import uuid4

from datetime import datetime

from monty.serialization import dumpfn, loadfn
from monty.dev import requires

from atomate.utils.utils import get_logger, env_chk
from atomate.vasp.analysis.lattice_dynamics import (
    fit_force_constants,
    get_cutoffs,
)
from atomate.vasp.analysis.phonopy import (
    get_phonon_dos,
    get_phonon_band_structure,
    get_line_mode_phonon_band_structure,
)
from atomate.vasp.database import VaspCalcDb
from atomate.vasp.workflows.base.lattice_dynamics import (
    MAX_N_IMAGINARY,
    MAX_IMAGINARY_FREQ,
    IMAGINARY_TOL,
    FIT_METHOD,
    DEFAULT_TEMPERATURE,
    MESH_DENSITY)
from fireworks import explicit_serialize, FiretaskBase, FWAction
import numpy as np

from pymatgen import Structure
from pymatgen.io.ase import AseAtomsAdaptor
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


@explicit_serialize
class CollectPerturbedStructures(FiretaskBase):
    """
    Aggregate the structures and forces of perturbed supercells. Requires a
    wf_uuid from which all perturbed static calculations will be compiled. The
    initial structure and supercell transformation matrix will be extracted from
    the perturbed structure tasks.

    Required parameters:
        db_file (str): Path to DB file for the database that contains the
            perturbed structure calculations.
        wf_uuid (str): Workflow unique identifier.
    """

    required_params = ["db_file", "wf_uuid"]

    def run_task(self, fw_spec):
        db_file = env_chk(self.get("db_file"), fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        tasks = mmdb.collection

        query = {
            "wf_meta.wf_uuid": self["wf_uuid"],
            "task_label": {"$regex": "perturbed static"},
            "state": "successful",
        }
        projection = ["output.structure", "output.forces", "transformations"]

        logger.info("Fetching structure and forces data")
        results = list(tasks.find(query, projection))

        logger.info("Found {} perturbed structures".format(len(results)))
        all_structures = []
        all_forces = []
        for result in results:
            structure = Structure.from_dict(result["output"]["structure"])
            forces = np.asarray(result["output"]["forces"])
            all_structures.append(structure)
            all_forces.append(forces)

        if len(all_structures) == 0:
            raise RuntimeError(
                "Could not find any perturbed structures for workflow uuid: "
                "{}".format(self["wf_uuid"])
            )
        else:
            transformation = results[-1]["transformations"][0]
            t_class = transformation["@class"]
            if t_class != "SupercellTransformation":
                raise ValueError(
                    "Expecting SupercellTransformation not {}".format(t_class)
                )
            structure = Structure.from_dict(transformation["input_structure"])
            supercell_matrix = transformation["deformation"]

        # regenerate pure supercell structure
        st = SupercellTransformation(supercell_matrix)
        supercell_structure = st.apply_transformation(structure)

        structure_data = {
            "structure": structure,
            "supercell_structure": supercell_structure,
            "supercell_matrix": supercell_matrix,
        }

        logger.info("Writing structure and forces data.")
        dumpfn(all_structures, "perturbed_structures.json")
        dumpfn(all_forces, "perturbed_forces.json")
        dumpfn(structure_data, "structure_data.json")


@requires(hiphive, "hiphive is required for lattice dynamics workflow")
@explicit_serialize
class FitForceConstants(FiretaskBase):
    """
    Used to aggregate atomic forces of perturbed supercells in and fit
    force constants using hiPhive. Requires "perturbed_structures.json",
    "perturbed_forces.json", and "structure_data.json" files to be present
    in the current working directory.

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

    def run_task(self, fw_spec):
        from hiphive.utilities import get_displacements

        max_n_imaginary = self.get("max_n_imaginary", MAX_N_IMAGINARY)
        max_imaginary_freq = self.get("max_imaginary_freq", MAX_IMAGINARY_FREQ)
        imaginary_tol = self.get("imaginary_tol", IMAGINARY_TOL)
        fit_method = self.get("fit_method", FIT_METHOD)

        all_structures = loadfn("perturbed_structures.json")
        all_forces = loadfn("all_forces.json")
        structure_data = loadfn("structure_data.json")
        parent_structure = structure_data["structure"]
        supercell_matrix = structure_data["supercell_matrix"]
        supercell_structure = structure_data["supercell_structure"]

        structures = []
        supercell_atoms = AseAtomsAdaptor.get_atoms(supercell_structure)
        for structure, forces in zip(all_structures, all_forces):
            atoms = AseAtomsAdaptor.get_atoms(structure)
            displacements = get_displacements(atoms, supercell_atoms)
            structure.new_array("displacements", displacements)
            structure.new_array("forces", forces)
            structure.positions = supercell_atoms.get_positions()
            structures.append(structure)

        cutoffs = self.get("cutoffs", get_cutoffs(supercell_structure))
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

        if force_constants is not None:
            # fitting failed
            raise RuntimeError(
                "Could not find a force constant solution with less than {}"
                "imaginary modes.\n"
                "Fitting results: {}".format(max_n_imaginary, fitting_data)
            )

        else:
            logger.info("Writing force constants.")
            force_constants.write("force_constants.fcs")

            atoms = AseAtomsAdaptor.get_atoms(parent_structure)
            force_constants.write_to_shengBTE("FORCE_CONSTANTS_3ND", atoms)
            force_constants.write_to_phonopy(
                "FORCE_CONSTANTS_2ND", format="text"
            )


@requires(hiphive, "hiphive is required for lattice dynamics workflow")
@explicit_serialize
class ForceConstantsToDB(FiretaskBase):
    """
    Add force constants, phonon band structure and density of states
    to the database.

    Assumes you are in a directory with the force constants, fitting
    data, and structure data written to files.

    Required parameters:
        db_file (str): Path to DB file for the database that contains the
            perturbed structure calculations.

    Optional parameters:
        wf_uuid (str): Workflow unique identifier.
        mesh_density (float): The density of the q-point mesh used to calculate
            the phonon density of states. See the docstring for the ``mesh``
            argument in Phonopy.init_mesh() for more details.
        additional_fields (dict): Additional fields added to the document, such
            as user-defined tags, name, ids, etc.
    """

    required_params = ["db_file"]
    optional_params = ["wf_uuid", "mesh_density", "additional_fields"]

    def run_task(self, fw_spec):
        from hiphive.force_constants import ForceConstants

        db_file = env_chk(self.get("db_file"), fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)

        force_constants = ForceConstants.read("force_constants.fcs")
        fitting_data = loadfn("fitting_data.json")
        structure_data = loadfn("structure_data.json")
        structure = structure_data["structure"]
        supercell_structure = structure_data["supercell_structure"]
        supercell_matrix = structure_data["supercell_matrix"]

        phonopy_fc = force_constants.get_fc_array(order=2)

        logger.info("Getting uniform phonon band structure.")
        uniform_bs = get_phonon_band_structure(
            structure, supercell_matrix, phonopy_fc
        )

        logger.info("Getting line mode phonon band structure.")
        lm_bs = get_line_mode_phonon_band_structure(
            structure, supercell_matrix, phonopy_fc
        )

        logger.info("Getting phonon density of states.")
        mesh_density = self.get("mesh_density", MESH_DENSITY)
        dos = get_phonon_dos(
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

        fc_uuid = str(uuid4())  # additional force constant identifier

        data = {
            "structure": structure.as_dict(),
            "supercell_matrix": supercell_matrix.tolist(),
            "supercell_structure": supercell_structure.as_dict(),
            "fitting_data": fitting_data,
            "force_constants": force_constants.get_fc_dict(),
            "tags": fw_spec.get("tags", None),
            "formula_pretty": structure.composition.reduced_formula,
            "phonon_dos_fs_id": dos_fsid,
            "phonon_bandstructure_uniform_fs_id": uniform_bs_fsid,
            "phonon_bandstructure_line_fs_id": lm_bs_fsid,
            "wf_meta": {"wf_uuid": self["wf_uuid"], "fc_uuid": fc_uuid},
            "created_at": datetime.utcnow(),
        }
        data.update(self.get("additional_fields", {}))

        mmdb.collection = mmdb.db["lattice_dynamics"]
        mmdb.collection.insert_one(data)

        return FWAction(update_spec={"fc_uuid": fc_uuid})


@explicit_serialize
class RunShengBTE(FiretaskBase):
    """
    Run ShengBTE to calculate lattice thermal conductivity. Presumes
    the FORCE_CONSTANTS_3ND and FORCE_CONSTANTS_2ND files in the current
    directory. In addition, presumes a "structure_data.json" file, with they
    keys "structure", " and "supercell_matrix" is in the current directory.

    Required parameters:
        shengbte_cmd (str): The name of the shengbte executable to run. Supports
            env_chk.

    Optional parameters:
        temperature (float or dict): The temperature to calculate the lattice
            thermal conductivity for. Can be given as a single float, or a
            dictionary with the keys "min", "max", "step".
        shengbte_control_kwargs (dict): Options to be included in the ShengBTE
            control file.
    """

    required_params = ["shengbte_cmd"]
    optional_params = ["temperature", "shengbte_control_kwargs"]

    def run_task(self, fw_spec):
        structure_data = loadfn("structure_data.json")
        structure = structure_data["structure"]
        supercell_matrix = structure_data["supercell_matrix"]

        control_dict = {
            "scalebroad": 0.5,
            "nonanalytic": False,
            "isotopes": False,
            "temperature": self.get("temperature", DEFAULT_TEMPERATURE),
            "scell": np.diag(supercell_matrix),
        }
        control_dict.update(self.get("shengbte_control_kwargs", {}))
        control = Control.from_structure(structure, **control_dict)
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
            return_code = os.subprocess.call(
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
class ShengBTEToDB(FiretaskBase):
    """
    Add lattice thermal conductivity results to database.

    Assumes you are in a directory with the ShengBTE results in.

    Required parameters:
        db_file (str): Path to DB file for the database that contains the
            perturbed structure calculations.

    Optional parameters:
        wf_uuid (str): Workflow unique identifier.
        fc_uuid (str): Force constant identifier, used if workflow has multiple
            force constants calculated (i.e., using different fit_methods).
        additional_fields (dict): Additional fields added to the document, such
            as user-defined tags, name, ids, etc.
    """

    required_params = ["db_file"]
    optional_params = ["wf_uuid", "fc_uuid", "additional_fields"]

    def run_task(self, fw_spec):
        db_file = env_chk(self.get("db_file"), fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)

        control = Control.from_file("CONTROL")
        structure = control.get_structure()
        supercell_matrix = np.diag(control["scell"])

        if Path("BTE.KappaTensorVsT_CONV").exists():
            filename = "BTE.KappaTensorVsT_CONV"
        elif Path("BTE.KappaTensorVsT_RTA").exists():
            filename = "BTE.KappaTensorVsT_RTA"
        else:
            raise RuntimeError("Could not find ShengBTE output files.")

        bte_data = np.loadtxt(filename)
        temperatures = bte_data[:, 0].tolist()
        kappa = bte_data[:, 1:10].reshape(-1, 3, 3)

        wf_meta = {}
        if self.get("wf_uuid"):
            wf_meta["wf_uuid"] = self["wf_uuid"]

        if self.get("fc_uuid") or fw_spec.get("fc_uuid"):
            wf_meta["fc_uuid"] = self.get("fc_uuid", fw_spec.get("fc_uuid"))

        data = {
            "structure": structure.as_dict(),
            "supercell_matrix": supercell_matrix.tolist(),
            "temperatures": temperatures,
            "lattice_thermal_conductivity": kappa,
            "control": control.to_json(),
            "tags": fw_spec.get("tags", None),
            "formula_pretty": structure.composition.reduced_formula,
            "wf_meta": wf_meta,
            "created_at": datetime.utcnow(),
        }
        data.update(self.get("additional_fields", {}))

        mmdb.collection = mmdb.db["lattice_thermal_conductivity"]
        mmdb.collection.insert(data)
