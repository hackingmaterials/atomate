from datetime import datetime
from json import loads

from fireworks import FiretaskBase, FWAction, explicit_serialize
from monty.json import MontyEncoder
from pydash import get
from pymatgen.analysis.path_finder import ChgcarPotential, NEBPathfinder
from pymatgen.core import Element, Structure
from pymatgen.io.vasp.sets import MPRelaxSet

from atomate.utils.utils import env_chk, get_logger
from atomate.vasp.database import VaspCalcDb

__author__ = "Ann Rutt"
__email__ = "acrutt@lbl.gov"

logger = get_logger(__name__)


@explicit_serialize
class HostToDb(FiretaskBase):
    """
    Initializes a approx_neb collection database entry from
    a host task doc specified by the supplied task_id and
    the provided approx_neb_wf_uuid. Also updates the
    tasks collection for approx neb workflow record keeping.

    Args:
        db_file (str): path to file containing the database
            credentials.
        approx_neb_wf_uuid (str): unique id for approx neb
            workflow record keeping
        host_task_id (int): task_id for VASP calculation of
            host structure. Must be provided in the fw_spec
            or firetask inputs.
    Optional Args:
        additional_fields (dict): specifies more information
            to be stored in the approx_neb collection to
            assist user record keeping.
        tags (list): list of strings to be stored in the
            approx_neb collection under the "tags" field to
            assist user record keeping.
    """

    required_params = ["db_file", "approx_neb_wf_uuid"]
    optional_params = ["host_task_id", "additional_fields", "tags"]

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
                "Provided approx_neb_wf_uuid is not unique. A unique workflow id is required "
                "for querying in the approx_neb workflow."
            )

        # update host task doc (from host_task_id) with unique wf_uuid
        # (tracks approx_neb workflows generated from this host task doc)
        t_id = self.get("host_task_id", fw_spec.get("host_task_id"))
        host_tasks_doc = mmdb.collection.find_one_and_update(
            {"task_id": t_id, "approx_neb.calc_type": "host"},
            {"$push": {"approx_neb.wf_uuids": wf_uuid}},
        )
        if host_tasks_doc is None:
            raise ValueError(f"Error updating approx neb host with task_id: {t_id}")

        # Initialize and store select host task doc fields in approx_neb_doc
        # (to be stored in the approx_neb collection)
        approx_neb_doc = {
            "wf_uuid": wf_uuid,
            "host": {
                "dir_name": host_tasks_doc["dir_name"],
                "chemsys": host_tasks_doc["chemsys"],
                "formula_pretty": host_tasks_doc["formula_pretty"],
                "input_structure": host_tasks_doc["input"]["structure"],
                "output": host_tasks_doc["output"],
                "task_id": host_tasks_doc["task_id"],
            },
            "end_points": [],
        }

        # ensure tags and additional_fields are the same
        # in both the approx_neb and tasks collections
        additional_fields = self.get("additional_fields", {})
        if isinstance(additional_fields, dict):
            for key, value in additional_fields.items():
                if key not in approx_neb_doc.keys():
                    approx_neb_doc[key] = value

        tags = self.get("tags")
        if tags:
            approx_neb_doc["tags"] = tags

        # insert approx_neb_doc in the approx_neb collection of provided database
        # includes fix to ensure approx_neb_doc is a json serializable dict
        approx_neb_doc = MontyEncoder().encode(approx_neb_doc)
        approx_neb_doc = loads(approx_neb_doc)
        approx_neb_doc["last_updated"] = datetime.utcnow()
        mmdb.collection = mmdb.db["approx_neb"]
        mmdb.collection.insert_one(approx_neb_doc)

        # Update fw_spec with approx_neb_doc and store wf_uuid
        # in launches collection for record keeping
        return FWAction(
            stored_data={"wf_uuid": wf_uuid, "approx_neb_doc": approx_neb_doc}
        )


@explicit_serialize
class PassFromDb(FiretaskBase):
    """
    Search approx_neb collection of database for designated
    fields. Pass designated fields to fw_spec according to
    fields_to_pull input.

    Args:
        db_file (str): path to file containing the database
            credentials.
        approx_neb_wf_uuid (str): unique identifier for
            approx neb workflow record keeping. Used for
            querying approx_neb collection.
        fields_to_pull (dict): define fields to pull from
            approx_neb collection using pydash.get()
            notation (doc specified by approx_neb_wf_uuid).
            Keys of fields_to_pull are used to name pulled
            fields in the updated fw_spec via...
            Format:
            {key : path} -> fw.spec[key] = task_doc[path]
            The path is a full mongo-style path so
            subdocuments can be referenced using dot
            notation and array keys can be referenced
            using the index. Example for pulling host
            structure from approx_neb collection into
            fw_spec["host_structure"]:
            {"host_structure":"host.output.structure"}
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
    Insert sites into the output structure of a previous
    calculation for use of the modified structure. Updates
    fw_spec with modified_structure to pass the modified
    structure. Also updates the fw_spec with
    end_points_index and structure_path if a
    approx_neb_wf_uuid is provided. Intended use is for
    inserting working ions into an empty host structure.

    Args:
        db_file (str): path to file containing the database
            credentials
        host_task_id (int): task_id for output structure to
            modify and insert site. Must be provided in the
            fw_spec or firetask inputs.
        insert_specie (str): specie of site to insert in
            structure (e.g. "Li").
        insert_coords (1x3 array or list of 1x3 arrays):
            coordinates of site(s) to insert in structure
            (e.g. [0,0,0] or [[0,0,0],[0,0.25,0]]).
        end_points_index (int): index used in end_points
            field of approx_neb collection for workflow
            record keeping.
        approx_neb_wf_uuid (str): Unique identifier for
            approx workflow record keeping.
    Optional Parameters:
        coords_are_cartesian (bool): Set to True if using
            cartesian coordinates for insert_coords.
            Otherwise assumes fractional coordinates.
    """

    required_params = [
        "db_file",
        "insert_specie",
        "insert_coords",
        "end_points_index",
        "approx_neb_wf_uuid",
    ]
    optional_params = ["host_task_id", "coords_are_cartesian"]

    def run_task(self, fw_spec):
        # get the database connection
        db_file = env_chk(self["db_file"], fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        wf_uuid = self.get("approx_neb_wf_uuid")

        insert_specie = self["insert_specie"]
        insert_coords = self["insert_coords"]
        end_points_index = self["end_points_index"]
        # put in list if insert_coords provided as a single coordinate to avoid error
        if isinstance(insert_coords[0], (float, int)):
            insert_coords = [insert_coords]

        # get output structure from host_task_id or approx_neb collection
        t_id = self.get("host_task_id", fw_spec.get("host_task_id"))
        if t_id:
            structure_doc = mmdb.collection.find_one({"task_id": t_id})
        else:
            mmdb.collection = mmdb.db["approx_neb"]
            approx_neb_doc = mmdb.collection.find_one({"wf_uuid": wf_uuid})
            structure_doc = approx_neb_doc["host"]["output"]["structure"]

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

        # store end point input structures in approx_neb collection
        mmdb.collection = mmdb.db["approx_neb"]
        mmdb.collection.update_one(
            {"wf_uuid": wf_uuid},
            {
                "$set": {
                    "end_points."
                    + str(end_points_index): {
                        "input_structure": structure.as_dict(),
                        "inserted_site_indexes": inserted_site_indexes,
                        "insert_coords": insert_coords,
                    }
                }
            },
        )

        stored_data = {
            "modified_structure": structure.as_dict(),
            "inserted_site_indexes": inserted_site_indexes,
            "structure_path": "end_points."
            + str(end_points_index)
            + ".input_structure",
        }

        return FWAction(
            update_spec={
                "insert_specie": insert_specie,
                "inserted_site_indexes": inserted_site_indexes,
            },
            stored_data=stored_data,
        )


@explicit_serialize
class WriteVaspInput(FiretaskBase):
    """
    Creates VASP input files using implementations of
    pymatgen's AbstractVaspInputSet. Vasp input parameters
    can be provided as a VaspInputSet object. The input
    structure used is set by the provided approx_neb_wf_uuid
    and structure_path to query the approx_neb collection
    using pydash.get().

    Args:
        db_file (str): Path to file specifying db
            credentials for getting the input structure from
            the approx_neb collection (specified by
            structure_path).
        approx_neb_wf_uuid (str): unique id for approx neb
            workflow record keeping.
        vasp_input_set (VaspInputSet class): can use to
            define VASP input parameters.
            See pymatgen.io.vasp.sets module for more
            information. MPRelaxSet() and
            override_default_vasp_params are used if
            vasp_input_set = None.
        structure_path (str): A full mongo-style path to
            reference approx_neb collection subdocuments
            using dot notation and array keys (e.g.
            "host.output.structure"). Must be provided in
            the fw_spec or firetask inputs.
        override_default_vasp_params (dict): if provided,
            vasp_input_set is disregarded and the Vasp Input
            Set is created by passing
            override_default_vasp_params to MPRelaxSet().
            Allows for easy modification of MPRelaxSet().
            For example, to set ISIF=2 in the INCAR use:
            {"user_incar_settings":{"ISIF":2}}
        potcar_spec (bool): Instead of writing the POTCAR, write a
            "POTCAR.spec". This is intended to allow testing of workflows
            without requiring pseudo-potentials to be installed on the system.
    """

    required_params = ["db_file", "approx_neb_wf_uuid", "vasp_input_set"]
    optional_params = [
        "structure_path",
        "override_default_vasp_params",
        "potcar_spec",
    ]

    def run_task(self, fw_spec):

        # get database connection
        db_file = env_chk(self["db_file"], fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        mmdb.collection = mmdb.db["approx_neb"]

        # get structure from approx_neb collection
        try:
            wf_uuid = self["approx_neb_wf_uuid"]
            structure_path = self.get("structure_path", fw_spec.get("structure_path"))
            approx_neb_doc = mmdb.collection.find_one({"wf_uuid": wf_uuid})
            structure = Structure.from_dict(get(approx_neb_doc, structure_path))

        except Exception:
            raise ValueError("Error getting structure from approx_neb collection")

        # get vasp input set and write files
        override_default_vasp_params = self.get("override_default_vasp_params") or {}
        if self["vasp_input_set"] is None or override_default_vasp_params != {}:
            vis = MPRelaxSet(
                structure, **override_default_vasp_params, sort_structure=False
            )
            # note sort_structure = False required to retain original site order of structure
        elif hasattr(self["vasp_input_set"], "write_input"):
            vis = self["vasp_input_set"]
        else:
            raise TypeError("ApproxNEB: Error using vasp_input_set")

        potcar_spec = self.get("potcar_spec", False)

        vis.write_input(".", potcar_spec=potcar_spec)


@explicit_serialize
class EndPointToDb(FiretaskBase):
    """
    Store information from VASP calculation of end point
    structure in approx_neb collection from the task doc
    specified by end_point_task_id. Also updates the tasks
    collection for approx neb workflow record keeping.

    Args:
        db_file (str): path to file containing the database
            credentials.
        approx_neb_wf_uuid (str): unique id for approx neb
            workflow record keeping
        end_points_index (int): index used in end_points
            field of approx_neb collection for workflow
            record keeping.
    Optional Params:
        end_point_task_id (int): task_id for VASP
            calculation of end point structure (empty host
            with one working ion inserted).
            end_point_task_id must be provided in the
            fw_spec or firetask inputs.
        wf_input_host_structure (structure dict): for record
            keeping purposes, a dictionary representation of
            the workflow's original input structure
        wf_input_insert_coords (list of 1x3 arrays): for
            record keeping purposes, the workflow's original
            input insert_coords (e.g. [[0,0,0]])
    """

    required_params = ["db_file", "approx_neb_wf_uuid", "end_points_index"]
    optional_params = [
        "end_point_task_id",
        "wf_input_host_structure",
        "wf_input_insert_coords",
    ]

    def run_task(self, fw_spec):
        # get the database connection
        db_file = env_chk(self["db_file"], fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        wf_uuid = self["approx_neb_wf_uuid"]
        index = self["end_points_index"]
        t_id = self.get(
            "end_point_task_id", fw_spec.get("end_points_" + str(index) + "_task_id")
        )
        print(t_id)

        # get any workflow inputs provided in fw_spec to store in task_doc
        wf_input_host_structure = fw_spec.get("wf_input_host_structure")
        wf_input_insert_coords = self.get("wf_input_insert_coords")
        wf_insertion_site_specie = fw_spec.get("insert_specie")
        wf_insertion_site_index = fw_spec.get("inserted_site_indexes")

        # get task doc (parts stored in approx_neb collection) and update for record keeping
        task_doc = mmdb.collection.find_one_and_update(
            {"task_id": t_id, "approx_neb.calc_type": "end_point"},
            {
                "$push": {
                    "approx_neb.wf_uuids": wf_uuid,
                    "approx_neb.end_points_indexes": index,
                },
                "$set": {
                    "approx_neb._wf_input_host_structure": wf_input_host_structure,
                    "approx_neb._wf_input_insert_coords": wf_input_insert_coords,
                    "approx_neb._wf_insertion_site_specie": wf_insertion_site_specie,
                    "approx_neb._wf_insertion_site_index": wf_insertion_site_index,
                },
            },
        )
        if task_doc is None:
            raise ValueError(
                f"Error updating approx neb end point with task_id: {t_id}"
            )

        # Store info in approx_neb collection for record keeping
        mmdb.collection = mmdb.db["approx_neb"]

        ep_subdoc = mmdb.collection.find_one(
            {"wf_uuid": wf_uuid},
            {"end_points": 1, "_id": 0},
        )
        ep_subdoc = ep_subdoc["end_points"]
        end_point_output = {
            "dir_name": task_doc["dir_name"],
            "formula_pretty": task_doc["formula_pretty"],
            "output": task_doc["output"],
            "task_id": task_doc["task_id"],
        }
        ep_subdoc[index].update(end_point_output)
        mmdb.collection.update_one(
            {"wf_uuid": wf_uuid},
            {"$set": {"end_points": ep_subdoc, "last_updated": datetime.utcnow()}},
        )

        return FWAction(
            stored_data={
                "wf_uuid": wf_uuid,
                "end_points_index": index,
                "end_point_output": end_point_output,
            }
        )


@explicit_serialize
class PathfinderToDb(FiretaskBase):
    """
    Applies NEBPathFinder (pymatgen.analysis.path_finder)
    using the host (CHGCAR from the task_id stored) and
    output structures stored in the end_points field of the
    approx_neb collection. Resulting image structures
    (interpolated between end point structures) are stored
    in the "images" field of the approx_neb collection for
    future use. The provided approx_neb_wf_uuid specifies
    the set of inputs to use.

    Args:
        db_file (str): path to file containing the database
            credentials.
        approx_neb_wf_uuid (str): unique id for approx neb
            workflow record keeping
        n_images: n_images (int): number of images
            interpolated between end point structures
    Optional Parameters:
        end_points_combo (str): string must have format of
            "0+1", "0+2", etc. to specify which combination
            of end_points to use for path interpolation.
    """

    required_params = ["db_file", "n_images", "approx_neb_wf_uuid", "end_points_combo"]
    optional_params = []

    def run_task(self, fw_spec):
        n_images = self["n_images"]
        end_points_combo = self["end_points_combo"]
        # checks if format of end_points_combo is correct and if so
        # get two desired end_points indexes for end point structures
        try:
            combo = end_points_combo.split("+")
            if len(combo) == 2:
                c = [int(combo[0]), int(combo[1])]
            else:
                raise ValueError("NEBPathfinder requires exactly two end points")
        except Exception:
            raise ValueError(
                f"{str(end_points_combo)} end_points_combo input is incorrect"
            )

        # get the database connection
        db_file = env_chk(self["db_file"], fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        mmdb.collection = mmdb.db["approx_neb"]
        wf_uuid = self["approx_neb_wf_uuid"]

        # get end_points and task_id of host from approx_neb collection
        approx_neb_doc = mmdb.collection.find_one(
            {"wf_uuid": wf_uuid}, {"end_points": 1, "host.task_id": 1, "_id": 0}
        )
        end_points = approx_neb_doc["end_points"]
        task_id = approx_neb_doc["host"]["task_id"]

        # get potential gradient, v, from host chgcar
        mmdb.collection = mmdb.db["tasks"]
        host_chgcar = mmdb.get_chgcar(task_id)
        v_chgcar = ChgcarPotential(host_chgcar)
        host_v = v_chgcar.get_v()

        # get start and end point structures from end_points
        start_struct = Structure.from_dict(end_points[c[0]]["output"]["structure"])
        end_struct = Structure.from_dict(end_points[c[1]]["output"]["structure"])

        # checks if inserted site indexes match
        inserted_site_indexes = end_points[c[0]]["inserted_site_indexes"]
        if inserted_site_indexes != end_points[c[1]]["inserted_site_indexes"]:
            raise ValueError(
                "Inserted site indexes of end point structures must match for NEBPathfinder"
            )

        # applies NEBPathFinder to interpolate and get images to store in
        # pathfinder_output.
        neb_pf = NEBPathfinder(
            start_struct,
            end_struct,
            relax_sites=inserted_site_indexes,
            v=host_v,
            n_images=n_images + 1,
        )
        # note NEBPathfinder currently returns n_images+1 images (rather than n_images)
        # and the first and last images generated are very similar to the end points
        # provided so they are discarded
        pathfinder_output = {
            "images": [structure.as_dict() for structure in neb_pf.images[1:-1]],
            "relax_site_indexes": inserted_site_indexes,
        }

        # stores images generated by NEBPathFinder in approx_neb collection
        # pathfinder field which is a nested dictionary using
        # end_points_combo as a key.
        mmdb.collection = mmdb.db["approx_neb"]
        pf_subdoc = mmdb.collection.find_one(
            {"wf_uuid": wf_uuid}, {"pathfinder": 1, "_id": 0}
        )
        if "pathfinder" not in pf_subdoc.keys():
            pf_subdoc = {}
        else:
            pf_subdoc = pf_subdoc["pathfinder"]
        pf_subdoc.update({end_points_combo: pathfinder_output})

        mmdb.collection.update_one(
            {"wf_uuid": wf_uuid},
            {"$set": {"pathfinder": pf_subdoc, "last_updated": datetime.utcnow()}},
        )

        return FWAction(
            stored_data={
                "wf_uuid": wf_uuid,
                "end_points_combo": c,
                "pathfinder": pathfinder_output,
            }
        )


@explicit_serialize
class AddSelectiveDynamics(FiretaskBase):
    """
    Applies specified selective dynamics scheme (see
    functions defined in Firetask) to interpolated images
    generated (and stored in approx_neb collection) by
    PathfinderToDb Firetask. Stores in the approx_neb
    collection (as input structures for future calculations).

    Args:
        db_file (str): path to file containing the database
            credentials
        approx_neb_wf_uuid (str): Unique identifier for
            approx workflow record keeping. Used to specify
            the approx_neb collection doc to get structures
            (from the "pathfinder" field), modify, and then
            store (in the "images" field).
        pathfinder_key (str): The "pathfinder" field of the
            approx_neb collection doc is a nested dictionary
            (to handle cases with multiple paths). This
            input specifies the key for this nested dict
            which is derived from the desired combination of
            end points. pathfinder_key should be a string
            of format "0+1", "0+2", etc. matching
            end_points_combo of the PathfinderToDb Firetask.
        mobile_specie (str): specie of site of interest such
            as the working ion (e.g. "Li" if the working ion
            of interest is lithium). Used to perform a check
            on the structures pulled from the approx_neb
            collection.
        selective_dynamics_scheme (str): "fix_two_atom"
    """

    required_params = [
        "db_file",
        "approx_neb_wf_uuid",
        "pathfinder_key",
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
            {"wf_uuid": wf_uuid}, {"pathfinder": 1, "images": 1}
        )

        # get structures stored in "pathfinder" field
        pathfinder_key = self["pathfinder_key"]
        structure_docs = approx_neb_doc["pathfinder"][pathfinder_key]["images"]
        fixed_index = approx_neb_doc["pathfinder"][pathfinder_key]["relax_site_indexes"]

        # apply selective dynamics to get images (list of pymatgen structures)
        fixed_specie = self["mobile_specie"]
        scheme = self["selective_dynamics_scheme"]
        images = self.get_images_list(structure_docs, scheme, fixed_index, fixed_specie)

        # assemble images output to be stored in approx_neb collection
        images_output = []
        for n, image in enumerate(images):
            images_output.append(
                {
                    "index": n,
                    "input_structure": image.as_dict(),
                    "selective_dynamics_scheme": scheme,
                }
            )

        # store images output in approx_neb collection
        if "images" not in approx_neb_doc.keys():
            images_subdoc = {}
        else:
            images_subdoc = approx_neb_doc["images"]
        images_subdoc.update({pathfinder_key: images_output})

        mmdb.collection.update_one(
            {"wf_uuid": wf_uuid},
            {"$set": {"images": images_subdoc, "last_updated": datetime.utcnow()}},
        )

        return FWAction(
            stored_data={"wf_uuid": wf_uuid, "images_output": images_output}
        )

    def get_images_list(self, structure_docs, scheme, fixed_index, fixed_specie):
        """
        Returns a list of structure objects with selective
        dynamics applied according to the specified scheme.

        Args:
            structure_docs (list): list of structure
                dictionaries
            scheme (str): "fix_two_atom"
        Returns:
             list of structure objects
        """
        if isinstance(structure_docs, (list)) is False:
            raise TypeError("list input required for structure_docs")

        if scheme == "fix_two_atom":
            images = []
            for doc in structure_docs:
                structure = Structure.from_dict(doc)
                image = self.add_fix_two_atom_selective_dynamics(
                    structure, fixed_index, fixed_specie
                )
                images.append(image)
        # ToDo: add radius based selective dynamics scheme
        else:
            raise ValueError(
                "selective_dynamics_scheme does match any supported schemes, check input value"
            )
        return images

    def add_fix_two_atom_selective_dynamics(self, structure, fixed_index, fixed_specie):
        """
        Returns structure with selective dynamics assigned
        to fix the position of two sites.
        Two sites will be fixed:
        1) the site specified by fixed_index and
        2) the site positioned furthest from the specified
           fixed_index site.

        Args:
            structure (Structure): Input structure (e.g.
                relaxed host with one working ion)
            fixed_index (int): Index of site in structure
                whose position will be fixed (e.g. the
                working ion site)
            fixed_specie (str or Element): Specie of site
                in structure whose position will be fixed
                (e.g. the working ion site)
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
    Store information from VASP calculation of image
    structure in approx_neb collection from the task doc
    specified by image_task_id. Also updates the tasks
    collection for approx neb workflow record keeping.

    Args:
        db_file (str): path to file containing the database
            credentials.
        approx_neb_wf_uuid (str): unique id for approx neb
            workflow record keeping
        image_task_id (int): task_id for VASP calculation
            of image structure (mobile specie/
            working ion in various positions in host).
            image_task_id must be provided in the fw_spec
            or firetask inputs.
        wf_input_host_structure (structure dict): for record
            keeping purposes, a dictionary representation of
            the workflow's original input structure
        wf_input_insert_coords (list of 1x3 arrays): for
            record keeping purposes, the workflow's original
            input insert_coords (e.g. [[0,0,0]])

        Note "approx_neb.image_index" field is required in
        task doc specified by image_task_id.
    """

    required_params = ["db_file", "approx_neb_wf_uuid"]
    optional_params = [
        "image_task_id",
        "wf_input_host_structure",
        "wf_input_insert_coords",
    ]

    def run_task(self, fw_spec):
        # get the database connection
        db_file = env_chk(self["db_file"], fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        wf_uuid = self["approx_neb_wf_uuid"]
        t_id = self.get("image_task_id", fw_spec.get("image_task_id"))

        # get any workflow inputs provided in fw_spec to store in task_doc
        wf_input_host_structure = fw_spec.get("wf_input_host_structure")
        wf_input_insert_coords = []
        for key, value in fw_spec.items():
            if "wf_input_insert_coords" in key:
                wf_input_insert_coords.append(value)

        # get task doc (parts stored in approx_neb collection) and update for record keeping
        task_doc = mmdb.collection.find_one_and_update(
            {"task_id": t_id, "approx_neb.calc_type": "image"},
            {
                "$push": {"approx_neb.wf_uuids": wf_uuid},
                "$set": {
                    "approx_neb._wf_input_host_structure": wf_input_host_structure,
                    "approx_neb._wf_input_insert_coords": wf_input_insert_coords,
                },
            },
        )
        if task_doc is None:
            raise ValueError(f"Error updating approx neb image with task_id: {t_id}")

        # Store info in approx_neb collection for record keeping
        images_key = task_doc["approx_neb"]["images_key"]
        index = task_doc["approx_neb"]["image_index"]
        mmdb.collection = mmdb.db["approx_neb"]
        images_subdoc = mmdb.collection.find_one(
            {"wf_uuid": wf_uuid}, {"images": 1, "_id": 0}
        )
        images_subdoc = images_subdoc["images"]
        image_output = {
            "dir_name": task_doc["dir_name"],
            "formula_pretty": task_doc["formula_pretty"],
            "output": task_doc["output"],
            "task_id": task_doc["task_id"],
        }
        images_subdoc[images_key][index].update(image_output)

        # path = "images." + images_key + "." + str(index) + "."
        # image_output = {
        #    path + "dir_name": task_doc["dir_name"],
        #    path + "formula_pretty": task_doc["formula_pretty"],
        #    path + "output": task_doc["output"],
        #    path + "task_id": task_doc["task_id"],
        #    "last_updated": datetime.utcnow(),
        # }

        mmdb.collection.update_one(
            {"wf_uuid": wf_uuid},
            {"$set": {"images": images_subdoc, "last_updated": datetime.utcnow()}},
        )

        return FWAction(
            stored_data={
                "wf_uuid": wf_uuid,
                "image_index": index,
                "image_output": image_output,
            }
        )
