from fireworks import Firework, FWAction, Workflow, FiretaskBase
from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.vasp.firetasks.glue_tasks import pass_vasp_result, CopyVaspOutputs
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet
from atomate.vasp.firetasks.parse_outputs import VaspToDb
from atomate.vasp.fireworks.core import OptimizeFW
from atomate.vasp.config import VASP_CMD, DB_FILE
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen import Structure, Element


# TODO: Write Firetask for taking CHGCAR from previous directory and storing it in database
@explicit_serialize
class HostLatticeToDb(FiretaskBase):
    """
    Initializes a approx_neb database entry from a host lattice task doc
    specified by the supplied task_id.
    Generates a unique id (wf_id) for approx_neb workflow record keeping.
    Host lattice task doc is updated with this unique approx_neb workflow id.
    Pulls GridFS identifier (or parses and stores from files in task doc
    directory if not already stored) for accessing host lattice CHGCAR and AECCARs.

    Required params:
        db_file (str): path to file containing the database credentials.
        host_lattice_task_id (int): task_id for structure optimization of host lattice
    """

    required_params = ["db_file", "host_lattice_task_id"]
    optional_params = []

    def run_task(self, fw_spec):
        # get the database connection
        db_file = env_chk(self.get("db_file"), fw_spec["db_file"])
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
        # and updates host lattice task doc with fs_id
        if any(fs_id == None for fs_id in [chgcar_fs_id, aeccar0_fs_id, aeccar2_fs_id]):
            calc_dir = approx_neb_doc["host_lattice"]["dir_name"]
            logger.info("CHECKING FOR CHG DENSITY FILES: {}".format(calc_dir))
            drone = VaspDrone(parse_chgcar=True, parse_aeccar=True)
            task_doc = drone.assimilate(calc_dir)
            # insert chgcar with GridFS
            if chgcar_fs_id == None:
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
                logger.info("CHGCAR GRIDFS INSERTION COMPLETE: {}".format(calc_dir))
            # insert aeccar with GridFS
            if aeccar0_fs_id == None or aeccar2_fs_id == None:
                aeccar0 = task_doc["calcs_reversed"][0]["aeccar0"]
                aeccar2 = task_doc["calcs_reversed"][0]["aeccar2"]
                # check if the aeccar is valid before insertion
                if (aeccar0.data["total"] + aeccar2.data["total"]).min() < 0:
                    logger.warning(
                        "AECCARs appear corrupted for task id {t_id} in {calc_dir}\nSkipping GridFS storage of AECCARs"
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
                    logger.info("AECCAR GRIDFS INSERTION COMPLETE: {}".format(calc_dir))

        # Store GridFS ids in approx_neb_doc (to be stored in the approx_neb collection)
        approx_neb_doc["chgcar_fs_id"] = chgcar_fs_id or chgcar_gfs_id
        approx_neb_doc["aeccar0_fs_id"] = aeccar0_fs_id or aeccar0_gfs_id
        approx_neb_doc["aeccar2_fs_id"] = aeccar2_fs_id or aeccar2_gfs_id
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
class InsertSites(FiretaskBase):
    """
    Insert sites into an output structure for use of the modified structure.

    Args:
        base_name (str): Name of the file to create.
        format (str): Optional. one of "zip", "tar", "bztar" or "gztar".
        """

    _fw_name = "InsertSites"
    required_params = ["insert_specie", "insert_coords"]
    optional_params = ["coords_are_cartesian"]

    def run_task(self, fw_spec):
        structure = Structure.from_dict(fw_spec["host_lattice_structure"])
        insert_coords = self["insert_coords"]
        insert_specie = self["insert_specie"]

        if structure.site_properties != {}:  # removes site properties to avoid error
            for p in structure.site_properties.keys():
                structure.remove_site_property(p)
        for coords in sorted(insert_coords, reverse=True):
            structure.insert(
                0, insert_specie, coords, **kwargs
            )  # add kwarg for coords_are_cartesian=False
        return structure


class InsertSitesFW(Firework):
    # TODO: Write class description
    def __init__(
        self,
        structure,
        insert_specie,
        insert_coords,
        name="approx neb insert working ion",
        vasp_input_set=None,
        override_default_vasp_params=None,
        vasp_cmd=VASP_CMD,
        db_file=DB_FILE,
        parents=None,
        **kwargs
    ):
        override_default_vasp_params = override_default_vasp_params or {}

        # if structure == None and parents == None:
        #   print("ERROR")
        # elif structure == None: #setting structure supercedes parents
        #   connect to database
        #   query for parent using fw_spec['_job_info'][-1]['launch_dir']
        #   get structure...
        #   structure #from parents
        # TODO:Is pass structure needed in this FW? How to ensure pass_dict key matches?
        pass_structure_fw = pass_vasp_result(
            pass_dict={"host_lattice_structure": ">>output.structure"}
        )
        structure = InsertSites(
            insert_specie=insert_specie, insert_coords=insert_coords
        )
        vasp_input_set = vasp_input_set or MPRelaxSet(
            structure, **override_default_vasp_params
        )
        t = [pass_structure_fw]
        t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, job_type="double_relaxation_run"))
        t.append(PassCalcLocs(name=name))
        t.append(VaspToDb(db_file=db_file, additional_fields={"task_label": name}))


class PathFinderFW(Firework):
    # TODO: Write PathfinderFW
    # FW requires starting from a previous calc to get CHGCAR
    def __init__(
        self,
        structure,
        insert_specie,
        insert_coords,
        parents=None,
        prev_calc_dir=None,
        name="pathfinder",
        vasp_input_set=None,
        override_default_vasp_params=None,
        vasp_cmd=VASP_CMD,
        db_file=DB_FILE,
        **kwargs
    ):
        override_default_vasp_params = override_default_vasp_params or {}
        t = []
        if prev_calc_dir:
            t.append(
                CopyVaspOutputs(calc_dir=prev_calc_dir, additional_files=["CHGCAR"])
            )
        elif parents:
            t.append(CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"]))
        else:
            raise ValueError(
                "Must specify previous calculation to use CHGCAR for PathfinderFW"
            )

        # ToDo: Apply Pathfinder
        task_name = name + "???"
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"))
        t.append(PassCalcLocs(name=task_name))
        t.append(VaspToDb(db_file=db_file, additional_fields={"task_label": task_name}))

    def add_fix_two_atom_selective_dynamics(structure, fixed_index, fixed_specie):
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
        if structure[fixed_index].specie != Element(fixed_specie):
            raise TypeError(
                "The chosen fixed atom at index {} is not a {} atom".format(
                    fixed_index, fixed_specie
                )
            )
        sd_structure = structure.copy()
        sd_array = [[True, True, True] for i in range(sd_structure.num_sites)]
        sd_array[fixed_index] = [False, False, False]
        ref_site = sd_structure.sites[fixed_index]
        distances = [site.distance(ref_site) for site in sd_structure.sites]
        farthest_index = distances.index(max(distances))
        sd_array[farthest_index] = [False, False, False]
        sd_structure.add_site_property("selective_dynamics", sd_array)
        return sd_structure
