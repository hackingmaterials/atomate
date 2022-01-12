from fireworks import FiretaskBase, FWAction, Workflow, explicit_serialize

from atomate.common.powerups import powerup_by_kwargs
from atomate.utils.utils import env_chk
from atomate.vasp.database import VaspCalcDb
from atomate.vasp.fireworks.approx_neb import ImageFW

__author__ = "Ann Rutt"
__email__ = "acrutt@lbl.gov"


@explicit_serialize
class GetImageFireworks(FiretaskBase):
    """
    Adds ImageFWs to the workflow for the provided images_key
    according to the scheme specified by launch_mode. Optional
    parameters such as "handler_group", "add_additional_fields",
    and "add_tags" can be used to modify the resulting ImageFWs.

    Args:
        db_file (str): path to file containing the database
            credentials.
        approx_neb_wf_uuid (str): unique id for approx neb workflow
            record keeping.
        images_key (str): specifies a key corresponding the images
            field of the approx_neb collection which specifies the
            desired combination of end points to interpolate images
            between. images_key should be a string of format "0+1",
            "0+2", etc. matching end_points_combo input of
            PathfinderToDb Firetask or pathfinder_key input of
            AddSelectiveDynamics Firetask. If images_key is not
            provided images will be launched for all paths/keys in
            the approx_neb collection images field.
        launch_mode (str): "all" or "screening"
        vasp_cmd (str): the name of the full executable for running
            VASP.
    Optional Params:
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
        handler_group (str or [ErrorHandler]): group of handlers to
            use for RunVaspCustodian firetask. See handler_groups
            dict in the code for the groups and complete list of
            handlers in each group. Alternatively, you can specify a
            list of ErrorHandler objects.
        add_additional_fields (dict): dict of additional fields to
            add to task docs (by additional_fields of VaspToDb).
        add_tags (list of strings): added to the "tags" field of the
            task docs.
    """

    required_params = [
        "db_file",
        "approx_neb_wf_uuid",
        "images_key",
        "launch_mode",
        "vasp_cmd",
    ]
    optional_params = [
        "vasp_input_set",
        "override_default_vasp_params",
        "handler_group",
        "add_additional_fields",
        "add_tags",
    ]

    def run_task(self, fw_spec):
        # get the database connection
        db_file = env_chk(self["db_file"], fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        mmdb.collection = mmdb.db["approx_neb"]
        wf_uuid = self["approx_neb_wf_uuid"]
        launch_mode = self["launch_mode"]
        images_key = self["images_key"]

        approx_neb_doc = mmdb.collection.find_one({"wf_uuid": wf_uuid}, {"images": 1})
        all_images = approx_neb_doc["images"]

        # get structure_path of desired images and sort into structure_paths
        if images_key and isinstance(all_images, (dict)):
            images = all_images[images_key]
            max_n = len(images)
            if launch_mode == "all":
                structure_paths = [
                    "images." + images_key + "." + str(n) + ".input_structure"
                    for n in range(0, max_n)
                ]
            elif launch_mode == "screening":
                structure_paths = self.get_and_sort_paths(
                    max_n=max_n, images_key=images_key
                )
        elif isinstance(all_images, (dict)):
            structure_paths = dict()
            if launch_mode == "all":
                for key, images in all_images.items():
                    max_n = len(images)
                    structure_paths[key] = [
                        "images." + key + "." + str(n) + ".input_structure"
                        for n in range(0, max_n)
                    ]
            elif launch_mode == "screening":
                for key, images in all_images.items():
                    structure_paths[key] = self.get_and_sort_paths(
                        max_n=len(images), images_key=key
                    )

        # get list of fireworks to launch
        if isinstance(structure_paths, list):
            if isinstance(structure_paths[0], (str)):
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
                    relax_image_fws.extend(
                        self.get_screening_fws(sorted_paths=sorted_paths)
                    )

        # place fws in temporary wf in order to use powerup_by_kwargs
        # to apply powerups to image fireworks
        if "vasp_powerups" in fw_spec.keys():
            temp_wf = Workflow(relax_image_fws)
            powerup_dicts = fw_spec["vasp_powerups"]
            temp_wf = powerup_by_kwargs(temp_wf, powerup_dicts)
            relax_image_fws = temp_wf.fws

        return FWAction(additions=relax_image_fws)

    def get_and_sort_paths(self, max_n, images_key=""):
        sorted_paths = [[], [], []]
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
        add_tags = self.get("add_tags")
        fw = ImageFW(
            approx_neb_wf_uuid=self["approx_neb_wf_uuid"],
            structure_path=structure_path,
            db_file=self["db_file"],
            vasp_input_set=self.get("vasp_input_set"),
            vasp_cmd=self["vasp_cmd"],
            override_default_vasp_params=self.get("override_default_vasp_params"),
            handler_group=self.get("handler_group"),
            parents=parents,
            add_additional_fields=self.get("add_additional_fields"),
            add_tags=add_tags,
        )
        if isinstance(add_tags, list):
            if "tags" in fw.spec.keys():
                fw.spec["tags"].extend(add_tags)
            else:
                fw.spec["tags"] = add_tags
        return fw

    def get_screening_fws(self, sorted_paths):
        if not isinstance(sorted_paths, list):
            if (
                any([isinstance(i, list) for i in sorted_paths]) is not True
                or len(sorted_paths) != 3
            ):
                raise TypeError("sorted_paths must be a list containing 3 lists")

        s1_fw = self.get_fw(structure_path=sorted_paths[0][0])
        # ToDo: modify this firework to add firetask that checks whether to run/defuse children

        s2_fws = []
        for path in sorted_paths[1]:
            s2_fws.append(self.get_fw(structure_path=path, parents=s1_fw))
            # ToDo: modify this firework to add firetask that checks whether to run/defuse children

        remaining_fws = []
        for path in sorted_paths[-1]:
            remaining_fws.append(self.get_fw(structure_path=path, parents=s2_fws))

        return [s1_fw] + s2_fws + remaining_fws
