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

@explicit_serialize
class GetImageFireworks(FiretaskBase):
    """
    ToDo: Update description

    Args:
        db_file (str): path to file containing the database credentials.
        approx_neb_wf_uuid (str): unique id for approx neb workflow record keeping
        launch_mode (int): "all" or "screening"
        vasp_cmd (...): ...
    Optional Params:
        images_key (str): for cases with multiple paths, to only launch images for
            one path use images_key to specify a key corresponding the images field
            derived from the desired combination of stable sites. images_key
            should be a string of format "0+1", "0+2", etc. matching
            stable_sites_combo input of PathfinderToDb Firetask or pathfinder_key
            input of AddSelectiveDynamics Firetask. If images_key is not provided
            images will be launched for all paths/keys in the approx_neb
            collection images field.
        handler_group (str or [ErrorHandler]): group of handlers to use for
            RunVaspCustodian firetask. See handler_groups dict in the code for
            the groups and complete list of handlers in each group. Alternatively,
            you can specify a list of ErrorHandler objects.
    """

    required_params = ["db_file", "approx_neb_wf_uuid", "launch_mode", "vasp_cmd"]
    optional_params = ["images_key", "vasp_input_set", "override_default_vasp_params", "handler_group", "add_additional_fields","add_tags"]

    def run_task(self, fw_spec):
        # get the database connection
        db_file = env_chk(self["db_file"], fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        mmdb.collection = mmdb.db["approx_neb"]
        wf_uuid = self["approx_neb_wf_uuid"]
        launch_mode = self["launch_mode"]
        images_key = self["images_key"]

        approx_neb_doc = mmdb.collection.find_one({"wf_uuid":wf_uuid},{"images":1})
        all_images = approx_neb_doc["images"]

        #get structure_path of desired images and sort into structure_paths
        if images_key and isinstance(all_images, (dict)):
            images = all_images[images_key]
            max_n = len(images)
            if launch_mode == "all":
                structure_paths = ["images." + images_key + "." + str(n) + ".input_structure" for n in range(0,max_n)]
            elif launch_mode == "screening":
                structure_paths = self.get_and_sort_paths(max_n=max_n, images_key=images_key)
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
            handler_group=self.get("handler_group"),
            parents = parents,
            add_additional_fields=self.get("add_additional_fields"),
            add_tags = self.get("add_tags")
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