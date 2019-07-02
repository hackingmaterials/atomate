from pymatgen.io.vasp.sets import MPRelaxSet
from fireworks import Firework
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet
from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.vasp.firetasks.parse_outputs import VaspToDb
from atomate.vasp.firetasks.approx_neb_tasks import (
    HostLatticeToDb,
    PassFromDb,
    InsertSites,
    WriteVaspInput,
    StableSiteToDb,
)
from atomate.vasp.config import VASP_CMD, DB_FILE


class HostLatticeFW(Firework):
    def __init__(
        self,
        structure,
        approx_neb_wf_uuid,
        name="approx neb host lattice",
        db_file=DB_FILE,
        vasp_input_set=None,
        vasp_cmd=VASP_CMD,
        override_default_vasp_params=None,
        job_type="double_relaxation_run",
        **kwargs
    ):
        """
        Launches a structure optimization calculation for a provided empty host lattice and
        stores appropriate fields in the task doc for approx_neb workflow record keeping.
        Stores initializes approx_neb collection database entry and stores relevant host
        lattice calcuation outputs.

        Adapted from OptimizeFW.

        Args:
            structure (Structure): input structure of empty host lattice
            approx_neb_wf_uuid (str): Unique identifier for approx workflow record keeping.
            name (str): name for the Firework.
            vasp_input_set (VaspInputSet): input set to use. Defaults to MPRelaxSet() if None.
            override_default_vasp_params (dict): If this is not None, these params are passed to
                the default vasp_input_set, i.e., MPRelaxSet. This allows one to easily override
                some settings, e.g., user_incar_settings, etc.
            vasp_cmd (str): Command to run vasp.
            db_file (str): Path to file specifying db credentials to place output parsing.
            job_type (str): custodian job type (default "double_relaxation_run")

            parents ([Firework]): Parents of this particular Firework.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        # set additional_fields to be added to task doc by VaspToDb
        # initiates the information stored in the tasks collection to aid record keeping
        additional_fields = {
            "task_label": name,
            "approx_neb": {
                "calc_type": "host_lattice",
                "wf_uuids": [],
                "_source_wf_uuid": approx_neb_wf_uuid,
            },
        }

        override_default_vasp_params = override_default_vasp_params or {}
        vasp_input_set = vasp_input_set or MPRelaxSet(
            structure, **override_default_vasp_params
        )

        t = []
        t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, job_type=job_type))
        t.append(PassCalcLocs(name=name))
        t.append(
            VaspToDb(
                db_file=db_file,
                additional_fields=additional_fields,
                parse_chgcar=True,
                parse_aeccar=True,
                task_fields_to_push={"host_lattice_task_id": "task_id"},
            )
        )
        t.append(
            HostLatticeToDb(db_file=db_file, approx_neb_wf_uuid=approx_neb_wf_uuid)
        )
        super().__init__(tasks=t, name="{} {}".format(structure.composition.reduced_formula, name), **kwargs)


class InsertSitesFW(Firework):
    def __init__(
        self,
        approx_neb_wf_uuid,
        insert_specie,
        insert_coords,
        name="approx neb stable site",
        db_file=DB_FILE,
        parents=None,
        **kwargs
    ):
        """
        Updates the fw_spec with the empty host lattice task_id from the provided
        approx_neb_wf_uuid. Pulls the empty host lattice structure from the tasks
        collection and inserts the site(s) designated by insert_specie and insert_coords.
        Stores the modified structure in the stable_sites field of the approx_neb collection.
        Updates the fw_spec with the corresponding stable_site_index for the stored
        structure (and the modified structure).

        Args:
            db_file (str): path to file containing the database credentials
            insert_specie (str): specie of site to insert in structure (e.g. "Li")
            insert_coords (1x3 array or list of 1x3 arrays): coordinates of site(s)
                to insert in structure (e.g. [0,0,0] or [[0,0,0],[0,0.25,0]])
            approx_neb_wf_uuid (str): Unique identifier for approx workflow record keeping.
            name (str): Name for the Firework.
            parents ([Firework]): Parents of this particular Firework.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        t = []
        # Add structure_task_id (for empty host lattice) to fw_spec
        # structure_task_id is required for InsertSites firetask
        t.append(
            PassFromDb(
                db_file=db_file,
                approx_neb_wf_uuid=approx_neb_wf_uuid,
                fields_to_pull={"structure_task_id": "host_lattice.task_id"},
            )
        )
        # Insert sites into empty host lattice (specified by structure_task_id)
        t.append(
            InsertSites(
                db_file=db_file,
                insert_specie=insert_specie,
                insert_coords=insert_coords,
                approx_neb_wf_uuid=approx_neb_wf_uuid,
            )
        )
        super().__init__(tasks=t, name=name, parents=parents, **kwargs)


class ApproxNEBLaunchFW(Firework):
    def __init__(
        self,
        calc_type,
        approx_neb_wf_uuid,
        structure_path=None,
        name="approx neb relaxation",
        db_file=DB_FILE,
        vasp_input_set=None,
        vasp_cmd=VASP_CMD,
        override_default_vasp_params=None,
        job_type="double_relaxation_run",
        parents=None,
        **kwargs
    ):
        """
        Launches a structure optimization calculation from a structure stored in
        in the approx_neb collection. Structure input for calculation is specified
        by the provided approx_neb_wf_uuid and structure_path to pull the structure
        from the approx_neb collection using pydash.get().

        Adapted from OptimizeFW.

        Args:
            calc_type (str): Set to
            approx_neb_wf_uuid (str): Unique identifier for approx workflow record keeping.
            structure_path (str): A full mongo-style path to reference approx_neb collection
                subdocuments using dot notation and array keys.
                By default structure_path = None which assumes fw_spec["structure_path"] is
                set by a parent firework.
            name (str): Name for the Firework.
            vasp_input_set (VaspInputSet): input set to use. Defaults to MPRelaxSet() if None.
            override_default_vasp_params (dict): If this is not None, these params are passed to
                the default vasp_input_set, i.e., MPRelaxSet. This allows one to easily override
                some settings, e.g., user_incar_settings, etc.
            vasp_cmd (str): Command to run vasp.
            db_file (str): Path to file specifying db credentials to place output parsing.
            job_type (str): custodian job type (default "double_relaxation_run")

            parents ([Firework]): Parents of this particular Firework.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """

        # set additional_fields to be added to task doc by VaspToDb
        # initiates the information stored in the tasks collection to aid record keeping
        additional_fields = {
            "task_label": name,
            "approx_neb": {"wf_uuids": [], "_source_wf_uuid": approx_neb_wf_uuid},
        }
        if calc_type == "stable_site":
            additional_fields["approx_neb"]["calc_type"] = "stable_site"
            additional_fields["approx_neb"]["stable_sites_indexes"]: []
        elif calc_type == "image":
            additional_fields["approx_neb"]["calc_type"] = "image"

        t = []
        t.append(
            WriteVaspInput(
                approx_neb_wf_uuid=approx_neb_wf_uuid,
                vasp_input_set=vasp_input_set,
                structure_path=structure_path,
                override_default_vasp_params=override_default_vasp_params,
            )
        )
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, job_type=job_type))
        t.append(PassCalcLocs(name=name))
        t.append(
            VaspToDb(
                db_file=db_file,
                additional_fields=additional_fields,
                task_fields_to_push={"stable_site_task_id": "task_id"},
            )
        )
        if calc_type == "stable_site":
            t.append(
                StableSiteToDb(db_file=db_file, approx_neb_wf_uuid=approx_neb_wf_uuid)
            )
        super().__init__(tasks=t, name=name, parents=parents, **kwargs)
