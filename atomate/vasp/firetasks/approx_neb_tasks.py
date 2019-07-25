from pymatgen import Element, Structure
from pymatgen.io.vasp.outputs import Chgcar
from pymatgen.analysis.path_finder import NEBPathfinder, ChgcarPotential
from fireworks import FiretaskBase, FWAction, explicit_serialize
from atomate.utils.utils import env_chk, load_class
from atomate.utils.utils import get_logger
from atomate.vasp.database import VaspCalcDb
from atomate.vasp.drones import VaspDrone
from pymatgen.io.vasp.sets import MPRelaxSet
from pydash import get  # ,set_
import json
import os.path
from monty.json import MontyEncoder
from itertools import combinations
from atomate.vasp.fireworks.approx_neb import ApproxNEBLaunchFW

logger = get_logger(__name__)


@explicit_serialize
class HostLatticeToDb(FiretaskBase):
    """
    Initializes a approx_neb collection database entry from a host lattice task doc
    specified by the supplied task_id using the provided approx_neb_wf_uuid.
    Host lattice task doc is updated with the provided approx_neb_wf_uuid.
    Pulls GridFS identifier (or parses and stores from files in task doc
    directory if not already stored) for accessing host lattice CHGCAR and AECCARs.

    Args:
        db_file (str): path to file containing the database credentials.
        approx_neb_wf_uuid (str): unique id for approx neb workflow record keeping
        host_lattice_task_id (int): task_id for structure optimization of host
            lattice. Must be provided in the fw_spec or firetask inputs.
    """

    required_params = ["db_file", "approx_neb_wf_uuid"]
    optional_params = ["host_lattice_task_id", "additional_fields"]

    def run_task(self, fw_spec):
        # get the database connection
        db_file = env_chk(self["db_file"], fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        wf_uuid = self["approx_neb_wf_uuid"]

        # check if provided approx_neb_wf_uuid is unique
        # e.g. not already used in approx_neb collection
        approx_neb_db = mmdb.db["approx_neb"]
        if approx_neb_db.count_documents({"wf_uuid": wf_uuid}) != 0:
            raise ValueError(
                "Provided approx_neb_wf_uuid is not unique. A unique workflow id is required for querying in the approx_neb workflow."
            )

        # get host lattice task doc from host_lattice_task_id
        t_id = self.get("host_lattice_task_id") or fw_spec.get("host_lattice_task_id")
        try:
            host_lattice_tasks_doc = mmdb.collection.find_one(
                {"task_id": t_id, "approx_neb.calc_type": "host_lattice"}
            )
        except:
            raise ValueError(
                "Unable to find approx neb host lattice with task_id: {}".format(t_id)
            )
        # update found host lattice task doc with unique wf_uuid
        # (tracks approx_neb workflows generated from this host lattice task doc)
        mmdb.collection.update_one(
            {"task_id": t_id}, {"$push": {"approx_neb.wf_uuids": wf_uuid}}
        )

        # Initialize and store select host lattice task doc fields in approx_neb_doc
        # (to be stored in the approx_neb collection)
        approx_neb_doc = {
            "host_lattice": {
                "dir_name": host_lattice_tasks_doc["dir_name"],
                "formula_pretty": host_lattice_tasks_doc["formula_pretty"],
                "output": host_lattice_tasks_doc["output"],
                "task_id": host_lattice_tasks_doc["task_id"],
            },
            "wf_uuid": wf_uuid,
            "stable_sites": [],
        }

        # Adds additional fields to approx_neb_doc stored in the task doc
        standard_task_doc_keys = [
            "_id",
            "dir_name",
            "task_label",
            "approx_neb",
            "schema",
            "calcs_reversed",
            "run_stats",
            "chemsys",
            "formula_anonymous",
            "formula_reduced_abc",
            "completed_at",
            "nsites",
            "composition_unit_cell",
            "composition_reduced",
            "formula_pretty",
            "elements",
            "nelements",
            "input",
            "output",
            "state",
            "analysis",
            "last_updated",
            "transformations",
            "custodian",
            "orig_inputs",
            "task_id",
        ]
        for key, value in host_lattice_tasks_doc.items():
            if key not in standard_task_doc_keys+["approx_neb"]:
                if key not in approx_neb_doc.keys():
                    approx_neb_doc[key] = value

        # Gets GridFS ids for host lattice chgcar and aeccar if stored in task_doc
        # gridfs_ids is for record keeping to track the source of GridFS ids
        chgcar_fs_id = host_lattice_tasks_doc["calcs_reversed"][0].get("chgcar_fs_id")
        aeccar0_fs_id = host_lattice_tasks_doc["calcs_reversed"][0].get("aeccar0_fs_id")
        aeccar2_fs_id = host_lattice_tasks_doc["calcs_reversed"][0].get("aeccar2_fs_id")
        gridfs_ids = {
            "task_doc_chgcar": chgcar_fs_id,
            "task_doc_aeccar0": aeccar0_fs_id,
            "task_doc_aeccar2": aeccar2_fs_id,
        }

        # If unable to get GridFS ids, checks task doc directory for
        # CHGCAR and/or AECCAR files. Stores files in database via GridFS
        # and updates host lattice task doc with the created fs_id.
        # Note: this will not work if computer does not have access to the
        # VASP calculation directory specified by the host lattice task doc dir_name
        # or if there is an error parsing the task doc dir_name
        if any(fs_id == None for fs_id in [chgcar_fs_id, aeccar0_fs_id, aeccar2_fs_id]):
            calc_dir = approx_neb_doc["host_lattice"]["dir_name"]
            # imperfect fix for parsing if host name is included task doc dir_name
            if ":" in calc_dir:
                calc_dir = calc_dir.split(":")[-1]

            logger.info("APPROX NEB: CHECKING FOR CHARGE DENSITY FILES")
            drone = VaspDrone(parse_chgcar=True, parse_aeccar=True)
            if os.path.exists(calc_dir):
                task_doc = drone.assimilate(calc_dir)
                output_files = task_doc["calcs_reversed"][0]["output_file_paths"].keys()

                # insert CHGCAR with GridFS
                if chgcar_fs_id == None and "chgcar" in output_files:
                    chgcar = json.dumps(
                        task_doc["calcs_reversed"][0]["chgcar"], cls=MontyEncoder
                    )
                    chgcar_gfs_id, compression_type = mmdb.insert_gridfs(
                        chgcar, "chgcar_fs", task_id=t_id
                    )
                    mmdb.collection.update_one(
                        {"task_id": t_id},
                        {
                            "$set": {
                                "calcs_reversed.0.chgcar_compression": compression_type,
                                "calcs_reversed.0.chgcar_fs_id": chgcar_gfs_id,
                            }
                        },
                    )
                    gridfs_ids["dir_chgcar"] = chgcar_gfs_id
                    logger.info("APPROX NEB: CHGCAR GRIDFS INSERTION COMPLETE")

                # insert AECCAR with GridFS
                if aeccar0_fs_id == None or aeccar2_fs_id == None:
                    if "aeccar0" in output_files and "aeccar2" in output_files:
                        aeccar0 = task_doc["calcs_reversed"][0]["aeccar0"]
                        aeccar2 = task_doc["calcs_reversed"][0]["aeccar2"]
                        # check if the aeccar is valid before insertion
                        if (aeccar0.data["total"] + aeccar2.data["total"]).min() < 0:
                            logger.warning(
                                "AECCARs APPEAR CORRUPTED. GridFS STORAGE SKIPPED"
                            )
                        else:
                            aeccar0_gfs_id, compression_type = mmdb.insert_gridfs(
                                aeccar0, "aeccar0_fs", task_id=t_id
                            )
                            mmdb.collection.update_one(
                                {"task_id": t_id},
                                {
                                    "$set": {
                                        "calcs_reversed.0.aeccar0_compression": compression_type,
                                        "calcs_reversed.0.aeccar0_fs_id": aeccar0_gfs_id,
                                    }
                                },
                            )
                            gridfs_ids["dir_aeccar0"] = aeccar0_gfs_id
                            aeccar2_gfs_id, compression_type = mmdb.insert_gridfs(
                                aeccar2, "aeccar2_fs", task_id=t_id
                            )
                            mmdb.collection.update_one(
                                {"task_id": t_id},
                                {
                                    "$set": {
                                        "calcs_reversed.0.aeccar2_compression": compression_type,
                                        "calcs_reversed.0.aeccar2_fs_id": aeccar2_gfs_id,
                                    }
                                },
                            )
                            gridfs_ids["dir_aeccar2"] = aeccar2_gfs_id
                            logger.info("APPROX NEB: AECCAR GRIDFS INSERTION COMPLETE")

        # Store GridFS ids in approx_neb_doc (to be stored in the approx_neb collection)
        # None will be stored if no gridfs_id is found
        approx_neb_doc["host_lattice"]["chgcar_fs_id"] = chgcar_fs_id or gridfs_ids.get("dir_chgcar")
        approx_neb_doc["host_lattice"]["aeccar0_fs_id"] = aeccar0_fs_id or gridfs_ids.get("dir_aeccar0")
        approx_neb_doc["host_lattice"]["aeccar2_fs_id"] = aeccar2_fs_id or gridfs_ids.get("dir_aeccar2")
        # Insert approx_neb_doc in the approx_neb collection of provided database
        mmdb.collection = mmdb.db["approx_neb"]
        mmdb.collection.insert_one(approx_neb_doc)

        # Update fw_spec with approx_neb_doc and
        # store wf_uuid and gridfs_ids in launches collection for record keeping
        return FWAction(
            stored_data={"wf_uuid": wf_uuid, "gridfs_ids": gridfs_ids},
            update_spec=approx_neb_doc,
        )


@explicit_serialize
class PassFromDb(FiretaskBase):
    """
        Search approx_neb collection of database for designated fields.
        Pass designated fields to fw_spec according to fields_to_pull input.

        Args:
            db_file (str): path to file containing the database credentials.
            approx_neb_wf_uuid (str): unique identifier for approx neb workflow
                record keeping. Used for querying approx_neb collection.
            fields_to_pull (dict): define fields to pull from approx_neb collection
                using pydash.get() notation (doc specified by approx_neb_wf_uuid).
                Keys of fields_to_pull are used to name pulled fields in the
                updated fw_spec via...
                Format: {key : path} -> fw.spec[key] = task_doc[path]
                The path is a full mongo-style path so subdocuments can be
                referenced using dot notation and array keys can be referenced
                using the index. Example for pulling host lattice structure from
                approx_neb collection into fw_spec["host_lattice_structure"]:
                {"host_lattice_structure":"host_lattice.output.structure"}
    """

    required_params = ["db_file", "approx_neb_wf_uuid", "fields_to_pull"]
    optional_params = []

    def run_task(self, fw_spec):
        # get the database connection
        db_file = env_chk(self["db_file"], fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        mmdb.collection = mmdb.db["approx_neb"]

        wf_uuid = self["approx_neb_wf_uuid"]
        fields_to_pull = self["fields_to_pull"]

        # pulls desired fields from approx_neb collection and stores in pulled_fields
        pulled_doc = mmdb.collection.find_one({"wf_uuid": wf_uuid})
        pulled_fields = dict()
        for key in fields_to_pull.keys():
            pulled_fields[key] = get(pulled_doc, fields_to_pull[key])

        # update fw_spec with pulled fields (labeled according to fields_to_pass)
        return FWAction(update_spec=pulled_fields)


@explicit_serialize
class InsertSites(FiretaskBase):
    """
    Insert sites into the output structure of a previous calculation
    for use of the modified structure. Updates fw_spec with
    modified_structure to pass the modified structure. Also updates the fw_spec
    with stable_sites_index and structure_path if a approx_neb_wf_uuid is provided.
    Intended use is for inserting working ions into an empty host lattice.

    Args:
        db_file (str): path to file containing the database credentials
        structure_task_id (int): task_id for output structure to modify and insert
            site. Must be provided in the fw_spec or firetask inputs.
        insert_specie (str): specie of site to insert in structure (e.g. "Li").
        insert_coords (1x3 array or list of 1x3 arrays): coordinates of site(s)
            to insert in structure (e.g. [0,0,0] or [[0,0,0],[0,0.25,0]]).
    Optional Parameters:
        approx_neb_wf_uuid (str): Unique identifier for approx workflow record
            keeping. If provided, checks if the output structure from
            structure_task_id matches the host lattice structure stored in the
            approx_neb collection doc specified by approx_neb_wf_uuid.
        coords_are_cartesian (bool): Set to True if you are providing insert_coords
            in cartesian coordinates. Otherwise assumes fractional coordinates.
        """

    required_params = ["db_file", "insert_specie", "insert_coords"]
    optional_params = [
        "structure_task_id",
        "approx_neb_wf_uuid",
        "coords_are_cartesian",
    ]

    def run_task(self, fw_spec):
        # get the database connection
        db_file = env_chk(self["db_file"], fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)

        insert_specie = self["insert_specie"]
        insert_coords = self["insert_coords"]
        # put in list if insert_coords provided as a single coordinate to avoid error
        if isinstance(insert_coords[0], (float, int)):
            insert_coords = [insert_coords]

        # get output structure from structure_task_id
        t_id = self.get("structure_task_id") or fw_spec.get("structure_task_id")
        try:
            structure_doc = mmdb.collection.find_one({"task_id": t_id})
        except:
            raise ValueError("Unable to find task doc with task_id: {}".format(t_id))

        # if provided, checks wf_uuid matches for approx_neb workflow record keeping
        if self.get("approx_neb_wf_uuid"):
            try:
                wf_uuid = self.get("approx_neb_wf_uuid")
                mmdb.collection = mmdb.db["approx_neb"]
                approx_neb_doc = mmdb.collection.find_one({"wf_uuid": wf_uuid})
                if (
                    structure_doc["output"]["structure"]
                    != approx_neb_doc["host_lattice"]["output"]["structure"]
                ):
                    raise ValueError(
                        "ApproxNEB: structure_task_id does not match approx_neb_wf_uuid host lattice"
                    )
            except:
                raise ValueError(
                    "ApproxNEB: Error matching structure_task_id and approx_neb_wf_uuid"
                )

        structure = Structure.from_dict(structure_doc["output"]["structure"])
        # removes site properties to avoid error
        if structure.site_properties != {}:
            for p in structure.site_properties.keys():
                structure.remove_site_property(p)

        # insert site(s) in structure and stores corresponding site index in inserted_site_indexes
        inserted_site_indexes = []
        for coords in sorted(insert_coords, reverse=True):
            structure.insert(
                0,
                insert_specie,
                coords,
                coords_are_cartesian=self.get("coords_are_cartesian", False),
            )
            inserted_site_indexes = [i + 1 for i in inserted_site_indexes]
            inserted_site_indexes.insert(0, 0)

        update_spec = {"modified_structure": structure.as_dict()}

        # store stable site input structures in approx_neb collection
        if self.get("approx_neb_wf_uuid"):
            try:
                # TODO: Store fw_id in approx_neb.stable_sites.input_structure
                mmdb.collection.update_one(
                    {"wf_uuid": wf_uuid},
                    {
                        "$push": {
                            "stable_sites": {
                                "input_structure": structure.as_dict(),
                                "inserted_site_indexes": inserted_site_indexes,
                            }
                        }
                    },
                )
                # get stable_sites_index to update the fw_spec for easier record keeping
                pulled = mmdb.collection.find_one(
                    {"wf_uuid": wf_uuid}, {"stable_sites"}
                )
                stable_sites_index = len(pulled["stable_sites"]) - 1
                update_spec["stable_sites_index"] = stable_sites_index
                update_spec["structure_path"] = (
                    "stable_sites." + str(stable_sites_index) + ".input_structure"
                )
            except:
                logger.warning("ApproxNEB: InsertSites FIRETASK STORING ERROR")

        return FWAction(update_spec=update_spec)


@explicit_serialize
class WriteVaspInput(FiretaskBase):
    """
    Creates VASP input files using implementations of pymatgen's
    AbstractVaspInputSet. Vasp input parameters can be provided as a VaspInputSet
    object. The input structure used is set by the provided approx_neb_wf_uuid and
    structure_path to query the approx_neb collection using pydash.get().

    Args:
        db_file (str): Path to file specifying db credentials for getting the input
            structure from the approx_neb collection (specified by structure_path).
        approx_neb_wf_uuid (str): unique id for approx neb workflow record keeping
        vasp_input_set (VaspInputSet class): can use to define VASP input
            parameters. See pymatgen.io.vasp.sets module for more information.
            MPRelaxSet() and override_default_vasp_params are used if
            vasp_input_set = None.
        structure_path (str): A full mongo-style path to reference approx_neb
            collection subdocuments using dot notation and array keys (e.g.
            "host_lattice.output.structure"). Must be provided in the fw_spec
            or firetask inputs.
        override_default_vasp_params (dict): if provided, vasp_input_set is
            disregarded and the Vasp Input Set is created by passing
            override_default_vasp_params to MPRelaxSet(). Allows for easy
            modification of MPRelaxSet(). For example, to set ISIF=2 in the INCAR
            use: override_default_vasp_params = {"user_incar_settings":{"ISIF":2}}
    """

    required_params = ["db_file", "approx_neb_wf_uuid", "vasp_input_set"]
    optional_params = ["structure_path", "override_default_vasp_params"]

    def run_task(self, fw_spec):

        # get database connection
        db_file = env_chk(self["db_file"], fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        mmdb.collection = mmdb.db["approx_neb"]

        # get structure from approx_neb collection
        try:
            wf_uuid = self["approx_neb_wf_uuid"]
            structure_path = self.get("structure_path") or fw_spec.get("structure_path")
            approx_neb_doc = mmdb.collection.find_one({"wf_uuid": wf_uuid})
            structure = Structure.from_dict(get(approx_neb_doc, structure_path))

        except:
            raise ValueError("Error getting structure from approx_neb collection")

        # get vasp input set and write files
        override_default_vasp_params = self.get("override_default_vasp_params") or {}
        if self["vasp_input_set"] == None or override_default_vasp_params != {}:
            vis = MPRelaxSet(structure, **override_default_vasp_params)
        elif hasattr(self["vasp_input_set"], "write_input"):
            vis = self["vasp_input_set"]
        else:
            raise TypeError("ApproxNEB: Error using vasp_input_set")
        vis.write_input(".")


@explicit_serialize
class StableSiteToDb(FiretaskBase):
    """
    Store information from stable site structure optimization in approx_neb
    collection database entry from a the task doc specified by stable_site_task_id.
    Also updates the tasks collection for approx neb workflow record keeping.

    Args:
        db_file (str): path to file containing the database credentials.
        approx_neb_wf_uuid (str): unique id for approx neb workflow record keeping
        stable_site_task_id (int): task_id for structure optimization of stable
            site structure (empty host lattice with one working ion inserted).
            stable_site_task_id must be provided in the fw_spec or firetask inputs.

        Note fw_spec["stable_sites_index"] is required.
    """

    required_params = ["db_file", "approx_neb_wf_uuid"]
    optional_params = ["stable_site_task_id"]

    def run_task(self, fw_spec):
        # get the database connection
        db_file = env_chk(self["db_file"], fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        wf_uuid = self["approx_neb_wf_uuid"]
        t_id = self.get("stable_site_task_id") or fw_spec.get("stable_site_task_id")

        # Note InsertSites firetask updates the the fw_spec with the stable_sites_index
        index = fw_spec["stable_sites_index"]

        # Store info in tasks collection for record keeping
        mmdb.collection.update_one(
            {"task_id": t_id, "approx_neb.calc_type": "stable_site"},
            {
                "$push": {
                    "approx_neb.wf_uuids": wf_uuid,
                    "approx_neb.stable_sites_indexes": index,
                }
            },
        )

        # pull task doc to store parts in approx_neb_collection
        task_doc = mmdb.collection.find_one(
            {"task_id": t_id, "approx_neb.calc_type": "stable_site"}
        )
        if task_doc == None:
            raise ValueError(
                "Unable to find approx neb stable site with task_id: {}".format(t_id)
            )

        # Store info in approx_neb collection for record keeping
        mmdb.collection = mmdb.db["approx_neb"]
        path = "stable_sites." + str(index) + "."
        stable_site_output = {
            path + "dir_name": task_doc["dir_name"],
            path + "formula_pretty": task_doc["formula_pretty"],
            path + "output": task_doc["output"],
            path + "task_id": task_doc["task_id"],
        }

        mmdb.collection.update_one({"wf_uuid": wf_uuid}, {"$set": stable_site_output})

        return FWAction(
            stored_data={
                "wf_uuid": wf_uuid,
                "stable_sites_index": index,
                "stable_site_output": stable_site_output,
            }
        )


@explicit_serialize
class PathfinderToDb(FiretaskBase):
    """
    Applies NEBPathFinder (from pymatgen.analysis.path_finder) using the
    host lattice (chgcar from the task_id stored) and output structures stored
    in the stable_sites field of the approx_neb collection. Resulting
    structures or images (interpolated between stable_site structures or end
    point structures) are stored in the "images" field of the approx_neb
    collection for future use. The provided approx_neb_wf_uuid specifies the
    set of inputs to use.

    Args:
        db_file (str): path to file containing the database credentials.
        approx_neb_wf_uuid (str): unique id for approx neb workflow record keeping
        n_images: n_images (int): number of images interpolated between end point structures
    """

    required_params = ["db_file", "n_images", "approx_neb_wf_uuid"]
    optional_params = []

    def run_task(self, fw_spec):
        n_images = self["n_images"]

        # get the database connection
        db_file = env_chk(self["db_file"], fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        mmdb.collection = mmdb.db["approx_neb"]
        wf_uuid = self["approx_neb_wf_uuid"]

        # get stable_sites and task_id of host lattice from approx_neb collection
        approx_neb_doc = mmdb.collection.find_one(
            {"wf_uuid": wf_uuid},
            {"stable_sites": 1, "host_lattice.task_id": 1, "_id": 0},
        )
        stable_sites = approx_neb_doc["stable_sites"]
        task_id = approx_neb_doc["host_lattice"]["task_id"]

        # get potential gradient, v, from host lattice chgcar
        mmdb.collection = mmdb.db["tasks"]
        host_lattice_chgcar = mmdb.get_chgcar(task_id)
        v_chgcar = ChgcarPotential(host_lattice_chgcar)
        host_lattice_v = v_chgcar.get_v()
        mmdb.collection = mmdb.db["approx_neb"]

        # get end point structures from stable_sites
        end_point_structs = [
            Structure.from_dict(stable_site["output"]["structure"])
            for stable_site in stable_sites
        ]

        # checks if all indexes match and if so, gets the inserted site indexes
        indexes = [i["inserted_site_indexes"] for i in approx_neb_doc["stable_sites"]]
        if all([i == indexes[0] for i in indexes]):  # if all indexes match
            inserted_site_indexes = indexes[0]
        else:
            ValueError(
                "Inserted site indexes of end point structures must match for NEBPathfinder"
            )

        # applies NEBPathFinder to interpolate and get images to store in
        # pathfinder_output. pathfinder_output is a unnested dictionary if
        # there are only two end point structures. pathfinder_output will be
        # a nested dictionary if there are more than two end point structures.
        n_end_point_structs = len(end_point_structs)
        if n_end_point_structs == 2:
            neb_pf = NEBPathfinder(
                end_point_structs[0],
                end_point_structs[-1],
                relax_sites=inserted_site_indexes,
                v=host_lattice_v,
                n_images=n_images - 1,
            )
            # note NEBPathfinder currently returns n_images+1 images (rather than n_images)

            pathfinder_output = {
                "images": [structure.as_dict() for structure in neb_pf.images],
                "relax_site_indexes": inserted_site_indexes,
            }
        elif n_end_point_structs > 2:
            pathfinder_output = dict()
            for c in combinations(range(0,n_end_point_structs),2):
                neb_pf = NEBPathfinder(
                    end_point_structs[c[0]],
                    end_point_structs[c[1]],
                    relax_sites=inserted_site_indexes,
                    v=host_lattice_v,
                    n_images=n_images - 1,
                )
                key = str(c[0]) + "+" + str(c[1])
                pathfinder_output[key] = {
                    "images": [structure.as_dict() for structure in neb_pf.images],
                    "relax_site_indexes": inserted_site_indexes,
                }
        else:
            ValueError(
                "NEBPathfinder requires at least two stable sites"
            )
            #"Firetask does not support NEBPathfinder with more than two stable sites"
            # ToDo: handling more than two provided stable sites

        # stores images generated by NEBPathFinder in approx_neb collection
        mmdb.collection.update_one(
            {"wf_uuid": wf_uuid}, {"$set": {"pathfinder": pathfinder_output}}
        )

        return FWAction(
            stored_data={
                "wf_uuid": wf_uuid,
                "num_end_point_structs": len(end_point_structs),
                "pathfinder": pathfinder_output,
            }
        )


@explicit_serialize
class AddSelectiveDynamics(FiretaskBase):
    """
    #ToDo: update description
    Applies specified selective dynamics schemes (see functions specified in this
    Firetask for descriptions) to image structures generated (and stored in
    approx_neb collection) by PathfinderToDb Firetask. Stores image structures in
    the approx_neb collection for as input structures for future relaxations.

    Args:
        db_file (str): path to file containing the database credentials
        approx_neb_wf_uuid (str): Unique identifier for approx workflow record
            keeping. Used to specify the approx_neb collection doc from which
            structure are pulled from the "pathfinder.images" field, modified, and
            then stored in the "images" field.
            images field.
        mobile_specie (str): specie of site of interest such as the working ion
            (e.g. "Li" if the working ion of interest is a Li). Provided  to
            perform a built in check on the structures pulled the approx_neb doc.
        selective_dynamics_scheme (str): "fix_two_atom"
        """

    required_params = [
        "db_file",
        "approx_neb_wf_uuid",
        "mobile_specie",
        "selective_dynamics_scheme",
    ]
    optional_params = []

    def run_task(self, fw_spec):
        # get the database connection to the approx_neb collection
        db_file = env_chk(self["db_file"], fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        mmdb.collection = mmdb.db["approx_neb"]

        # get approx_neb doc specified by wf_uuid
        wf_uuid = self.get("approx_neb_wf_uuid")
        approx_neb_doc = mmdb.collection.find_one(
            {"wf_uuid": wf_uuid}, {"pathfinder": 1}
        )

        fixed_specie = self["mobile_specie"]
        scheme = self["selective_dynamics_scheme"]
        # get structures stored in "pathfinder" field
        # apply selective dynamics to get images
        if "images" in approx_neb_doc["pathfinder"].keys():
            structure_docs = approx_neb_doc["pathfinder"]["images"]
            fixed_index = approx_neb_doc["pathfinder"]["relax_site_indexes"]
            all_images = self.get_images_list(structure_docs,scheme,fixed_index,fixed_specie)
        else:
            all_images = dict()
            for key in approx_neb_doc["pathfinder"].keys():
                structure_docs = approx_neb_doc["pathfinder"][key]["images"]
                fixed_index = approx_neb_doc["pathfinder"][key]["relax_site_indexes"]
                images = self.get_images_list(structure_docs,scheme,fixed_index,fixed_specie)
                all_images[key] = images

        # assemble images output to be stored in approx_neb collection
        if isinstance(all_images, (list)):
            images_output = []
            for n, image in enumerate(all_images):
                images_output.append(
                    {
                        "index": n,
                        "input_structure": image.as_dict(),
                        "selective_dynamics_scheme": scheme,
                    }
                )
        elif isinstance(all_images, (dict)):
            images_output = dict()
            for key,images in all_images.items():
                images_output[key] = []
                for n, image in enumerate(images):
                    images_output[key].append({
                        "index": n,
                        "input_structure": image.as_dict(),
                        "selective_dynamics_scheme": scheme,
                    })

        mmdb.collection.update_one(
            {"wf_uuid": wf_uuid}, {"$set": {"images": images_output}}
        )

        return FWAction(
            update_spec={"images": images_output},
            stored_data={"wf_uuid": wf_uuid, "images_output": images_output},
        )

    def get_images_list(self, structure_docs, scheme, fixed_index, fixed_specie):
        """
        Returns a list of structure objects with selective dynamics applied
        according ot the specified scheme.

        Args:
            structure_docs(list of structure dictionaries):
            scheme(str): "fix_two_atom" or
        Returns:
             list of structure objects
        """
        if isinstance(structure_docs, (list)) == False:
            raise TypeError("list input required for structure_docs")

        if scheme == "fix_two_atom":
            images = []
            for doc in structure_docs:
                structure = Structure.from_dict(doc)
                image = self.add_fix_two_atom_selective_dynamics(
                    structure, fixed_index, fixed_specie
                )
                images.append(image)
        else:
            raise ValueError(
                "selective_dynamics_scheme does match any supported schemes, check input value"
            )
            # ToDo: add radius based selective dynamics scheme
        return images

    def add_fix_two_atom_selective_dynamics(self, structure, fixed_index, fixed_specie):
        """
        Returns structure with selective dynamics assigned to fix the
        position of two sites.
        Two sites will be fixed: 1) the site specified by fixed_index and
        2) the site positioned furthest from the specified fixed_index site.

        Args:
            structure (Structure): Input structure (e.g. host lattice with
            one working ion intercalated)
            fixed_index (int): Index of site in structure whose position
            will be fixed (e.g. the working ion site)
            fixed_specie (str or Element): Specie of site in structure
            whose position will be fixed (e.g. the working ion site)
        Returns:
            Structure
        """

        if isinstance(fixed_index, (list)):
            # avoids error if single index is provided as a list
            if len(fixed_index) == 1:
                fixed_index = fixed_index[0]
            else:
                raise ValueError(
                    "fixed_index must specify exactly one index for the fix_two_atom selective dynamics scheme"
                )

        if structure[fixed_index].specie != Element(fixed_specie):
            raise ValueError(
                "The chosen fixed atom at index {} is not a {} atom".format(
                    fixed_index, fixed_specie
                )
            )

        # removes site properties to avoid error
        if structure.site_properties != {}:
            for p in structure.site_properties.keys():
                structure.remove_site_property(p)

        sd_structure = structure.copy()
        sd_array = [[True, True, True] for i in range(sd_structure.num_sites)]
        sd_array[fixed_index] = [False, False, False]
        ref_site = sd_structure.sites[fixed_index]
        distances = [site.distance(ref_site) for site in sd_structure.sites]
        farthest_index = distances.index(max(distances))
        sd_array[farthest_index] = [False, False, False]
        sd_structure.add_site_property("selective_dynamics", sd_array)

        return sd_structure


@explicit_serialize
class ImageToDb(FiretaskBase):
    """
    ToDo: Update description
    Store information from stable site structure optimization in approx_neb
    collection database entry from a the task doc specified by stable_site_task_id.
    Also updates the tasks collection for approx neb workflow record keeping.

    Args:
        db_file (str): path to file containing the database credentials.
        approx_neb_wf_uuid (str): unique id for approx neb workflow record keeping
        image_task_id (int): task_id for structure optimization of image structure
            (mobile specie/working ion in various positions moving in host lattice
            structure).
            image_task_id must be provided in the fw_spec or firetask inputs.

        Note task_doc["approx_neb.image_index"] is required for task doc specified
        by image_task_id.
    """

    required_params = ["db_file", "approx_neb_wf_uuid"]
    optional_params = ["image_task_id"]

    def run_task(self, fw_spec):
        # get the database connection
        db_file = env_chk(self["db_file"], fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        wf_uuid = self["approx_neb_wf_uuid"]
        t_id = self.get("image_task_id") or fw_spec.get("image_task_id")

        # Store info in tasks collection for record keeping
        mmdb.collection.update_one(
            {"task_id": t_id, "approx_neb.calc_type": "image"},
            {"$push": {"approx_neb.wf_uuids": wf_uuid}},
        )

        # pull task doc to store parts in approx_neb_collection
        task_doc = mmdb.collection.find_one(
            {"task_id": t_id, "approx_neb.calc_type": "image"}
        )
        if task_doc == None:
            raise ValueError(
                "Unable to find approx neb image with task_id: {}".format(t_id)
            )

        # Store info in approx_neb collection for record keeping
        index = task_doc["approx_neb"]["image_index"]
        mmdb.collection = mmdb.db["approx_neb"]
        if "images_key" in task_doc["approx_neb"].keys():
            images_key = task_doc["approx_neb"]["images_key"]
            path = "images." + images_key + "." + str(index) + "."
        else:
            path = "images." + str(index) + "."
        image_output = {
            path + "dir_name": task_doc["dir_name"],
            path + "formula_pretty": task_doc["formula_pretty"],
            path + "output": task_doc["output"],
            path + "task_id": task_doc["task_id"],
        }

        mmdb.collection.update_one({"wf_uuid": wf_uuid}, {"$set": image_output})

        return FWAction(
            stored_data={
                "wf_uuid": wf_uuid,
                "image_index": index,
                "image_output": image_output,
            }
        )


@explicit_serialize
class LaunchImages(FiretaskBase):
    """
    ToDo: Update description

    Args:
        db_file (str): path to file containing the database credentials.
        approx_neb_wf_uuid (str): unique id for approx neb workflow record keeping
        launch_mode (int): "all" or "screening"
    """

    required_params = ["db_file", "approx_neb_wf_uuid", "launch_mode", "vasp_cmd"]
    optional_params = ["vasp_input_set", "override_default_vasp_params"]

    def run_task(self, fw_spec):
        # get the database connection
        db_file = env_chk(self["db_file"], fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        mmdb.collection = mmdb.db["approx_neb"]
        wf_uuid = self["approx_neb_wf_uuid"]
        launch_mode = self["launch_mode"]

        approx_neb_doc = mmdb.collection.find_one({"wf_uuid":wf_uuid},{"images":1})
        all_images = approx_neb_doc["images"]

        #get structure_path of desired images and sort into structure_paths
        if isinstance(all_images, (list)):
            max_n = len(all_images)
            if launch_mode == "all":
                structure_paths = ["images." + str(n) + ".input_structure" for n in range(0,max_n)]
            elif launch_mode == "screening":
                structure_paths = self.get_and_sort_paths(max_n)
        elif isinstance(all_images, (dict)):
            structure_paths = dict()
            if launch_mode == "all":
                for key, images in all_images.items():
                    max_n = len(images)
                    structure_paths[key] = ["images." + key + "." + str(n) + ".input_structure" for n in range(0,max_n)]
            elif launch_mode == "screening":
                for key, images in all_images.items():
                    structure_paths[key] = self.get_and_sort_paths(max_n=len(images), images_key=key)

        # get list of fireworks to launch
        if isinstance(structure_paths, (list)):
            if isinstance(structure_paths[0],(str)):
                relax_image_fws = []
                for path in structure_paths:
                    relax_image_fws.append(self.get_fw(structure_path=path))
            else:
                relax_image_fws = self.get_screening_fws(sorted_paths=structure_paths)
        elif isinstance(structure_paths, (dict)):
            relax_image_fws = []
            if launch_mode == "all":
                for key in structure_paths.keys():
                    for path in structure_paths[key]:
                        relax_image_fws.append(self.get_fw(structure_path=path))
            elif launch_mode == "screening":
                for key in structure_paths.keys():
                    sorted_paths = structure_paths[key]
                    relax_image_fws.extend(self.get_screening_fws(sorted_paths=sorted_paths))

        return FWAction(additions = relax_image_fws)

    def get_and_sort_paths(self, max_n, images_key=""):
        sorted_paths = [[],[],[]]

        mid_n = int(max_n / 2)
        q1 = int((max_n - mid_n) / 2)  # for second screening pass
        q3 = int((max_n + mid_n) / 2)  # for second screening pass

        for n in range(0, max_n):
            path = "images." + images_key + "." + str(n) + ".input_structure"
            if n == mid_n:  # path for first screening pass (center image index)
                sorted_paths[0].append(path)
            elif n in [q1, q3]:
                sorted_paths[1].append(path)
            else:
                sorted_paths[-1].append(path)

        return sorted_paths

    def get_fw(self, structure_path, parents=None):
        fw = ApproxNEBLaunchFW(
            calc_type="image",
            approx_neb_wf_uuid=self["approx_neb_wf_uuid"],
            structure_path=structure_path,
            db_file=self["db_file"],
            vasp_input_set=self.get("vasp_input_set"),
            vasp_cmd=self["vasp_cmd"],
            override_default_vasp_params=self.get("override_default_vasp_params"),
            parents = parents
        )
        return fw

    def get_screening_fws(self, sorted_paths):
        if isinstance(sorted_paths,(list))!=True:
            if any([isinstance(i,(list)) for i in sorted_paths]) != True or len(sorted_paths)!=3:
                raise TypeError("sorted_paths must be a list containing 3 lists")

        s1_fw = self.get_fw(structure_path=sorted_paths[0][0])
        # ToDo: modify this firework to add firetask that checks whether to run/defuse children

        s2_fws = []
        for path in sorted_paths[1]:
            s2_fws.append(self.get_fw(structure_path=path, parents=s1_fw))
            #ToDo: modify this firework to add firetask that checks whether to run/defuse children

        remaining_fws = []
        for path in sorted_paths[-1]:
            remaining_fws.append(self.get_fw(structure_path=path,parents=s2_fws))

        return [s1_fw] + s2_fws + remaining_fws















