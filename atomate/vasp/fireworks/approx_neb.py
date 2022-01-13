from fireworks import Firework
from pymatgen.io.vasp.sets import MPRelaxSet

from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.vasp.config import DB_FILE, VASP_CMD
from atomate.vasp.firetasks.approx_neb_tasks import (
    EndPointToDb,
    HostToDb,
    ImageToDb,
    InsertSites,
    PassFromDb,
    WriteVaspInput,
)
from atomate.vasp.firetasks.parse_outputs import VaspToDb
from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet

__author__ = "Ann Rutt"
__email__ = "acrutt@lbl.gov"


class HostFW(Firework):
    def __init__(
        self,
        structure,
        approx_neb_wf_uuid,
        db_file=DB_FILE,
        vasp_input_set=None,
        vasp_cmd=VASP_CMD,
        override_default_vasp_params=None,
        job_type="double_relaxation_run",
        additional_fields=None,
        tags=None,
        **kwargs,
    ):
        r"""
        Launches a VASP calculation for the provided empty host structure
        and adds task doc fields for approx_neb workflow record keeping.
        Initializes a doc in the approx_neb collection and stores relevant
        outputs from the host.

        Args:
            structure (Structure): structure of empty host
            approx_neb_wf_uuid (str): unique id for approx neb
                workflow record keeping
            db_file (str): path to file containing the database
                credentials.
            vasp_input_set (VaspInputSet class): can use to
                define VASP input parameters.
                See pymatgen.io.vasp.sets module for more
                information. MPRelaxSet() and
                override_default_vasp_params are used if
                vasp_input_set = None.
            vasp_cmd (str): the name of the full executable for running
                VASP.
            override_default_vasp_params (dict): if provided,
                vasp_input_set is disregarded and the Vasp Input
                Set is created by passing
                override_default_vasp_params to MPRelaxSet().
                Allows for easy modification of MPRelaxSet().
                For example, to set ISIF=2 in the INCAR use:
                {"user_incar_settings":{"ISIF":2}}
            job_type (str): custodian job type
            additional_fields (dict): specifies more information
                to be stored in the approx_neb collection to
                assist user record keeping.
            tags (list): list of strings to be stored in the
                approx_neb collection under the "tags" field to
                assist user record keeping.

            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
            parents ([Firework]): Parents of this particular Firework.
        """
        # set additional_fields to be added to task doc by VaspToDb
        # initiates the information stored in the tasks collection to aid record keeping
        fw_name = "{} {}".format(structure.composition.reduced_formula, "host")
        fw_spec = {"tags": ["approx_neb", approx_neb_wf_uuid, "host", "relaxation"]}
        task_doc_additional_fields = {
            "approx_neb": {
                "wf_uuids": [],
                "_source_wf_uuid": approx_neb_wf_uuid,
                "calc_type": "host",
                "task_label": "relaxation",
            }
        }

        override_default_vasp_params = override_default_vasp_params or {}
        vasp_input_set = vasp_input_set or MPRelaxSet(
            structure, **override_default_vasp_params
        )

        t = []
        t.append(WriteVaspFromIOSet(structure=structure, vasp_input_set=vasp_input_set))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, job_type=job_type))
        t.append(PassCalcLocs(name="host"))
        t.append(
            VaspToDb(
                db_file=db_file,
                additional_fields=task_doc_additional_fields,
                parse_chgcar=True,
                task_fields_to_push={"host_task_id": "task_id"},
            )
        )
        t.append(
            HostToDb(
                db_file=db_file,
                approx_neb_wf_uuid=approx_neb_wf_uuid,
                additional_fields=additional_fields,
                tags=tags,
            )
        )
        super().__init__(tasks=t, spec=fw_spec, name=fw_name, **kwargs)


class EndPointFW(Firework):
    def __init__(
        self,
        approx_neb_wf_uuid,
        insert_specie,
        insert_coords,
        end_points_index,
        db_file=DB_FILE,
        vasp_input_set=None,
        vasp_cmd=VASP_CMD,
        override_default_vasp_params=None,
        job_type="double_relaxation_run",
        parents=None,
        **kwargs,
    ):
        r"""
        Pulls information from the approx_neb collection (e.g. host
        task_id) using the provided approx_neb_wf_uuid. Gets the host
        structure from the tasks collection and inserts the site(s)
        designated by insert_specie and insert_coords. Stores the modified
        structure in the end_points field of the approx_neb collection.
        Launches a VASP calculation and adds task doc fields for approx_neb
        workflow record keeping. The input structure is specified by the
        provided approx_neb_wf_uuid and end_points_index. Stores relevant
        outputs in the approx_neb collection.

        Args:
            approx_neb_wf_uuid (str): unique id for approx neb
                workflow record keeping
            insert_specie (str): specie of site to insert in
                structure (e.g. "Li").
            insert_coords (1x3 array or list of 1x3 arrays):
                fractional coordinates of site(s) to insert in
                structure (e.g. [0,0,0] or [[0,0,0],[0,0.25,0]]).
            end_points_index (int): index used in end_points
                field of approx_neb collection for workflow
                record keeping.
            db_file (str): path to file containing the database
                credentials
            vasp_input_set (VaspInputSet class): can use to
                define VASP input parameters.
                See pymatgen.io.vasp.sets module for more
                information. MPRelaxSet() and
                override_default_vasp_params are used if
                vasp_input_set = None.
            vasp_cmd (str): the name of the full executable for running
                VASP.
            override_default_vasp_params (dict): if provided,
                vasp_input_set is disregarded and the Vasp Input
                Set is created by passing
                override_default_vasp_params to MPRelaxSet().
                Allows for easy modification of MPRelaxSet().
                For example, to set ISIF=2 in the INCAR use:
                {"user_incar_settings":{"ISIF":2}}
            job_type (str): custodian job type
            parents ([Firework]): Parents of this particular Firework.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        fw_name = f"end point: insert {insert_specie} {end_points_index}"
        fw_spec = {
            "tags": ["approx_neb", approx_neb_wf_uuid, "end_point", "relaxation"]
        }

        # set additional_fields to be added to task doc by VaspToDb
        # initiates the information stored in the tasks collection to aid record keeping
        additional_fields = {
            "approx_neb": {
                "wf_uuids": [],
                "_source_wf_uuid": approx_neb_wf_uuid,
                "_wf_input_host_structure": None,
                "_wf_input_insert_coords": None,
                "_wf_insertion_site_specie": None,
                "_wf_insertion_site_index": None,
                "calc_type": "end_point",
                "task_label": "relaxation",
                "end_points_indexes": [],
            }
        }

        t = []
        # Add host_task_id (for empty host) to fw_spec
        # host_task_id is required for InsertSites firetask.
        # wf_input_host_structure and wf_input_insert_coords are pulled
        # into the fw_spec for EndPointToDb firetask.
        t.append(
            PassFromDb(
                db_file=db_file,
                approx_neb_wf_uuid=approx_neb_wf_uuid,
                fields_to_pull={
                    "host_task_id": "host.task_id",
                    "wf_input_host_structure": "host.input_structure",
                    "wf_input_insert_coords_"
                    + str(end_points_index): "end_points."
                    + str(end_points_index)
                    + ".insert_coords",
                },
            )
        )
        # Insert sites into empty host (specified by host_task_id)
        t.append(
            InsertSites(
                db_file=db_file,
                insert_specie=insert_specie,
                insert_coords=insert_coords,
                end_points_index=end_points_index,
                approx_neb_wf_uuid=approx_neb_wf_uuid,
            )
        )
        # write vasp inputs, run vasp, parse vasp outputs
        t.append(
            WriteVaspInput(
                db_file=db_file,
                approx_neb_wf_uuid=approx_neb_wf_uuid,
                vasp_input_set=vasp_input_set,
                structure_path="end_points."
                + str(end_points_index)
                + ".input_structure",
                override_default_vasp_params=override_default_vasp_params,
            )
        )
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd, job_type=job_type))
        t.append(PassCalcLocs(name="end_point"))
        t.append(
            VaspToDb(
                db_file=db_file,
                additional_fields=additional_fields,
                task_fields_to_push={
                    "end_points_" + str(end_points_index) + "_task_id": "task_id"
                },
            )
        )
        # store desired outputs from tasks doc in approx_neb collection
        t.append(
            EndPointToDb(
                end_points_index=end_points_index,
                db_file=db_file,
                approx_neb_wf_uuid=approx_neb_wf_uuid,
                wf_input_insert_coords=insert_coords,
            )
        )

        super().__init__(tasks=t, spec=fw_spec, name=fw_name, parents=parents, **kwargs)


class ImageFW(Firework):
    def __init__(
        self,
        approx_neb_wf_uuid,
        structure_path=None,
        db_file=DB_FILE,
        vasp_input_set=None,
        vasp_cmd=VASP_CMD,
        override_default_vasp_params=None,
        job_type="double_relaxation_run",
        handler_group=None,
        parents=None,
        add_additional_fields=None,
        add_tags=None,
        **kwargs,
    ):
        r"""
        Pulls information from the approx_neb collection using the
        provided approx_neb_wf_uuid (including the image structure using
        structure_path and pydash.get() notation). Launches a VASP
        calculation and adds task doc fields for approx_neb workflow
        record keeping. Stores relevant outputs in the approx_neb
        collection.

        Args:
            approx_neb_wf_uuid (str): Unique identifier for approx_neb
                workflow record keeping.
            structure_path (str): A full mongo-style path to reference
                approx_neb collection subdocuments using dot notation and
                array keys. e.g. "images.0+1.2.input_structure"
                By default structure_path = None which assumes
                fw_spec["structure_path"] is set by a parent firework.
            vasp_input_set (VaspInputSet class): can use to
                define VASP input parameters.
                See pymatgen.io.vasp.sets module for more
                information. MPRelaxSet() and
                override_default_vasp_params are used if
                vasp_input_set = None.
            override_default_vasp_params (dict): if provided,
                vasp_input_set is disregarded and the Vasp Input
                Set is created by passing
                override_default_vasp_params to MPRelaxSet().
                Allows for easy modification of MPRelaxSet().
                For example, to set ISIF=2 in the INCAR use:
                {"user_incar_settings":{"ISIF":2}}
            vasp_cmd (str): the name of the full executable for running
                VASP.
            db_file (str): Path to file specifying db credentials to store
                outputs.
            job_type (str): custodian job type
            handler_group (str or [ErrorHandler]): group of handlers to use
                for RunVaspCustodian firetask. See handler_groups dict in
                the code for the groups and complete list of handlers in
                each group. Alternatively, you can specify a list of
                ErrorHandler objects.
            parents ([Firework]): Parents of this particular Firework.
            add_additional_fields (dict): dict of additional fields to
                add to task docs (by additional_fields of VaspToDb).
            add_tags (list of strings): added to the "tags" field of the
                task docs.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        # initiates the information stored in the tasks collection to aid record keeping
        image_index = int(structure_path.split(".")[-2])
        images_key = structure_path.split(".")[-3]
        fw_name = "image " + images_key + ": " + str(image_index)
        fw_spec = {"tags": ["approx_neb", approx_neb_wf_uuid, "image", "relaxation"]}
        handler_group = handler_group or {}

        # set additional_fields to be added to task doc by VaspToDb
        additional_fields = {
            "approx_neb": {
                "wf_uuids": [],
                "_source_wf_uuid": approx_neb_wf_uuid,
                "_wf_input_host_structure": None,
                "_wf_input_insert_coords": None,
                "calc_type": "image",
                "task_label": "relaxation",
                "images_key": images_key,
                "image_index": image_index,
            }
        }
        if isinstance(add_additional_fields, (dict)):
            for key, value in add_additional_fields.items():
                additional_fields[key] = value
        if isinstance(add_tags, (list, dict)):
            additional_fields["tags"] = add_tags

        t = []
        # pull input host and insert_coords for ImageToDb firetask
        combo = images_key.split("+")
        if len(combo) == 2:
            c = [int(combo[0]), int(combo[1])]
        t.append(
            PassFromDb(
                db_file=db_file,
                approx_neb_wf_uuid=approx_neb_wf_uuid,
                fields_to_pull={
                    "wf_input_host_structure": "host.input_structure",
                    "wf_input_insert_coords"
                    + str(c[0]): "end_points."
                    + str(c[0])
                    + ".insert_coords",
                    "wf_input_insert_coords"
                    + str(c[1]): "end_points."
                    + str(c[1])
                    + ".insert_coords",
                },
            )
        )
        # write vasp inputs, run vasp, parse vasp outputs
        t.append(
            WriteVaspInput(
                db_file=db_file,
                approx_neb_wf_uuid=approx_neb_wf_uuid,
                vasp_input_set=vasp_input_set,
                structure_path=structure_path,
                override_default_vasp_params=override_default_vasp_params,
            )
        )
        t.append(
            RunVaspCustodian(
                vasp_cmd=vasp_cmd, job_type=job_type, handler_group=handler_group
            )
        )
        t.append(PassCalcLocs(name="image"))
        t.append(
            VaspToDb(
                db_file=db_file,
                additional_fields=additional_fields,
                task_fields_to_push={"image_task_id": "task_id"},
            )
        )
        # store desired outputs from tasks doc in approx_neb collection
        t.append(ImageToDb(db_file=db_file, approx_neb_wf_uuid=approx_neb_wf_uuid))

        super().__init__(tasks=t, spec=fw_spec, name=fw_name, parents=parents, **kwargs)
