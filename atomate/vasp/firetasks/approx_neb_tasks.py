from pymatgen import Structure
from fireworks import FiretaskBase, FWAction, explicit_serialize
from atomate.utils.utils import env_chk
from atomate.utils.utils import get_logger
from atomate.vasp.database import VaspCalcDb
from atomate.vasp.drones import VaspDrone
from pydash import get  # ,set_
import json
import os.path
from monty.json import MontyEncoder
from uuid import uuid4

logger = get_logger(__name__)


@explicit_serialize
class HostLatticeToDb(FiretaskBase):
    """
    Initializes a approx_neb database entry from a host lattice task doc
    specified by the supplied task_id.
    Generates a unique id (wf_id) for approx_neb workflow record keeping.
    Host lattice task doc is updated with this unique approx_neb workflow id.
    Pulls GridFS identifier (or parses and stores from files in task doc
    directory if not already stored) for accessing host lattice CHGCAR and AECCARs.

    Args:
        db_file (str): path to file containing the database credentials.
        host_lattice_task_id (int): task_id for structure optimization of host lattice
    """

    required_params = ["db_file", "host_lattice_task_id"]
    optional_params = []

    def run_task(self, fw_spec):
        # get the database connection
        db_file = env_chk(self["db_file"], fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        wf_uuid = str(uuid4())

        # get host lattice task doc from host_lattice_task_id
        t_id = self["host_lattice_task_id"]
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
        }

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
        approx_neb_doc["chgcar_fs_id"] = chgcar_fs_id or gridfs_ids.get("dir_chgcar")
        approx_neb_doc["aeccar0_fs_id"] = aeccar0_fs_id or gridfs_ids.get("dir_aeccar0")
        approx_neb_doc["aeccar2_fs_id"] = aeccar2_fs_id or gridfs_ids.get("dir_aeccar2")
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
            approx_neb_wf_uuid (str): unique id for approx neb workflow record keeping
            fields_to_pull (dict): define fields to pull from approx_neb collection
             using pydash.get() notation (doc specified by approx_neb_wf_uuid).
                Keys of fields_to_pull are used to name pulled fields in the
                updated fw_spec. Example input (pull output structure into fw_spec):
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
    "modified_structure" field to pass the modified structure.
    Intended use is for inserting working ions into an empty host lattice.

    Args:
        db_file (str): path to file containing the database credentials
        structure_task_id (int): task_id for output structure to modify and insert site
        insert_specie (str): specie of site to insert in structure (e.g. "Li")
        insert_coords (1x3 array or list of 1x3 arrays): coordinates of site(s)
            to insert in structure (e.g. [0,0,0] or [[0,0,0],[0,0.25,0]])
    Optional Parameters:
        coords_are_cartesian (bool): Set to True if you are providing insert_coords
            in cartesian coordinates. Otherwise assumes fractional coordinates.
        approx_neb_wf_uuid (str): Unique identifier for approx workflow record keeping.
            If provided, checks if the output structure from structure_task_id matches
            the host lattice structure stored in the approx_neb collection doc specified
            by approx_neb_wf_uuid.
        """

    required_params = ["db_file", "structure_task_id", "insert_specie", "insert_coords"]
    optional_params = ["coords_are_cartesian", "approx_neb_wf_uuid"]

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
        t_id = self["structure_task_id"]
        try:
            structure_doc = mmdb.collection.find_one({"task_id": t_id})
        except:
            raise ValueError("Unable to find tasks doc with task_id: {}".format(t_id))

        # if provided, checks wf_uuid matches for approx_neb workflow record keeping
        if self.get("approx_neb_wf_uuid"):
            try:
                wf_uuid = self.get("approx_neb_wf_uuid")
                mmdb.collection = mmdb.db["approx_neb"]
                approx_neb_doc = mmdb.collection.find_one({"wf_uuid": wf_uuid})
                if structure_doc["output"]["structure"] != approx_neb_doc["host_lattice"]["output"]["structure"]:
                    raise ValueError("ApproxNEB: structure_task_id does not match approx_neb_wf_uuid host lattice")
            except:
                raise ValueError("ApproxNEB: Error matching structure_task_id and approx_neb_wf_uuid")

        structure = Structure.from_dict(structure_doc["output"]["structure"])
        # removes site properties to avoid error
        if structure.site_properties != {}:
            for p in structure.site_properties.keys():
                structure.remove_site_property(p)

        # insert site(s) in structure
        for coords in sorted(insert_coords, reverse=True):
            structure.insert(
                0,
                insert_specie,
                coords,
                coords_are_cartesian=self.get("coords_are_cartesian", False),
            )

        return FWAction(update_spec={"modified_structure": structure.to_json()})
