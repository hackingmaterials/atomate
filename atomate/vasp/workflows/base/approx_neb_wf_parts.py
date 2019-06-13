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


#TODO: Write Firetask for taking CHGCAR from previous directory and storing it in database
@explicit_serialize
class ChgDensityToDb(FiretaskBase):
    """
    Enter the CHGCAR or AECCAR from a previous run into a database.
    Uses current directory unless you specify the previous calculation using calc_dir or calc_loc.
    Adapted from VaspToDb Firetask but only stores charge density data.

    Required params:
        db_file (str): path to file containing the database credentials.
        parse_chgcar (bool):
        parse_aeccar (bool):
    Optional params:
        calc_dir (str): path to dir (on current filesystem) that contains VASP
            output files. Default: use current working directory.
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
            search for the most recent calc_loc with the matching name
        additional_fields (dict): dict of additional fields to add
        fw_spec_field (str): if set, will update the task doc with the contents
            of this key in the fw_spec.
        defuse_unsuccessful (bool): this is a three-way toggle on what to do if
            your job looks OK, but is actually unconverged (either electronic or
            ionic). True -> mark job as COMPLETED, but defuse children.
            False --> do nothing, continue with workflow as normal. "fizzle"
            --> throw an error (mark this job as FIZZLED)
        task_fields_to_push (dict): if set, will update the next Firework/Firetask
            spec using fields from the task document.
            Format: {key : path} -> fw.spec[key] = task_doc[path]
            The path is a full mongo-style path so subdocuments can be referneced
            using dot notation and array keys can be referenced using the index.
            E.g "calcs_reversed.0.output.outar.run_stats"
    """
    required_params = ["db_file", "parse_chgcar", "parse_aeccar"]
    optional_params = ["calc_dir", "calc_loc",
                       "additional_fields", "fw_spec_field", "defuse_unsuccessful",
                       "task_fields_to_push"]

    def run_task(self, fw_spec):
        # get the directory that contains the previous VASP calculation for storing CHGCAR and/or AECCAR
        calc_dir = os.getcwd()
        if "calc_dir" in self:
            calc_dir = self["calc_dir"]
        elif self.get("calc_loc"):
            calc_dir = get_calc_loc(self["calc_loc"], fw_spec["calc_locs"])["path"]

        # parse the VASP directory
        logger.info("PARSING DIRECTORY: {}".format(calc_dir))

        drone = VaspDrone(additional_fields=self.get("additional_fields"),
                          parse_chgcar=fw_spec["parse_chgcar"],
                          parse_aeccar=fw_spec["parse_aeccar"])

        # assimilate (i.e. parse)
        task_doc = drone.assimilate(calc_dir)

        # Check for additional keys to set based on the fw_spec
        if self.get("fw_spec_field"):
            task_doc.update(fw_spec[self.get("fw_spec_field")])

        # get the database connection
        db_file = fw_spec["db_file"] #env_chk(self.get('db_file'), fw_spec)

        # db insertion
        if db_file:
            mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
            #the chgcar
            if fw_spec["parse_chgcar"]:
                task_doc["calcs_reversed"][0]["chgcar"]

            if fw_spec["parse_aeccar"]:
                task_doc["calcs_reversed"][0]["aeccar0"]
                task_doc["calcs_reversed"][0]["aeccar2"]

            t_id = mmdb.insert_task(
                task_doc, self.get("parse_chgcar", False)
                or self.get("parse_aeccar", False))
            logger.info("Finished parsing with task_id: {}".format(t_id))

        defuse_children = False
        if task_doc["state"] != "successful":
            defuse_unsuccessful = self.get("defuse_unsuccessful",
                                           DEFUSE_UNSUCCESSFUL)
            if defuse_unsuccessful is True:
                defuse_children = True
            elif defuse_unsuccessful is False:
                pass
            elif defuse_unsuccessful == "fizzle":
                raise RuntimeError(
                    "ChgDensityToDb indicates that job is not successful "
                    "(perhaps your job did not converge within the "
                    "limit of electronic/ionic iterations)!")
            else:
                raise RuntimeError("Unknown option for defuse_unsuccessful: "
                                   "{}".format(defuse_unsuccessful))

        task_fields_to_push = self.get("task_fields_to_push", None)
        update_spec = {}
        if task_fields_to_push:
            if isinstance(task_fields_to_push, dict):
                for key, path_in_task_doc in task_fields_to_push.items():
                    if has(task_doc, path_in_task_doc):
                        update_spec[key] = get(task_doc, path_in_task_doc)
                    else:
                        logger.warning("Could not find {} in task document. Unable to push to next firetask/firework".format(path_in_task_doc))
            else:
                raise RuntimeError("Inappropriate type {} for task_fields_to_push. It must be a "
                                   "dictionary of format: {key: path} where key refers to a field "
                                   "in the spec and path is a full mongo-style path to a "
                                   "field in the task document".format(type(task_fields_to_push)))

        return FWAction(stored_data={"task_id": task_doc.get("task_id", None)},
                        defuse_children=defuse_children, update_spec=update_spec)

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

        if (
            structure.site_properties != {}
        ):  # removes site properties to avoid error
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
        t.append(
            WriteVaspFromIOSet(
                structure=structure, vasp_input_set=vasp_input_set
            )
        )
        t.append(
            RunVaspCustodian(
                vasp_cmd=vasp_cmd, job_type="double_relaxation_run"
            )
        )
        t.append(PassCalcLocs(name=name))
        t.append(
            VaspToDb(db_file=db_file, additional_fields={"task_label": name})
        )


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
                CopyVaspOutputs(
                    calc_dir=prev_calc_dir, additional_files=["CHGCAR"]
                )
            )
        elif parents:
            t.append(
                CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"])
            )
        else:
            raise ValueError(
                "Must specify previous calculation to use CHGCAR for PathfinderFW"
            )

        # ToDo: Apply Pathfinder
        task_name = name + "???"
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, auto_npar=">>auto_npar<<"))
        t.append(PassCalcLocs(name=task_name))
        t.append(
            VaspToDb(
                db_file=db_file, additional_fields={"task_label": task_name}
            )
        )

    def add_fix_two_atom_selective_dynamics(
        structure, fixed_index, fixed_specie
    ):
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
