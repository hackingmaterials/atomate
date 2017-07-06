# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import glob
import shutil

from pymatgen.core import Structure
from pymatgen.io.vasp import Incar, Kpoints, Poscar, Potcar
from pymatgen_diffusion.neb.io import MVLCINEBSet, get_endpoint_dist, get_endpoints_from_index

from fireworks.core.firework import FiretaskBase, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize

from atomate.utils.utils import get_logger

"""
NEB workflow firetasks.
"""

__author__ = "Hanmei Tang, Iek-Heng Chu"
__email__ = 'hat003@eng.ucsd.edu, ihchu@eng.ucsd.edu'


VASP_NEB_OUTPUT_FILES = {'INCAR', 'KPOINTS', 'POTCAR', 'vasprun.xml'}
VASP_NEB_OUTPUT_SUB_FILES = {'CHG', 'CHGCAR', 'CONTCAR', 'DOSCAR', 'EIGENVAL', 'IBZKPT', 'PCDAT',
                             'POSCAR', 'PROCAR', 'OSZICAR', 'OUTCAR', 'REPORT', 'WAVECAR', 'XDATCAR'}


@explicit_serialize
class TransferNEBTask(FiretaskBase):
    """
    This class transfers NEB outputs from current directory to destination directory. "label" is
    used to determine the step of calculation and hence the final path. The corresponding structure
    will be updated in fw_spec before files transferring.

    Required params:
        label (str): Type of calculation outputs being transferred, choose from "parent", "ep0",
            "ep1", "neb1", "neb2" and etc..
    Optional params:
        d_img (float): Distance between neighbouring images, used to determine the number of images
            if "IMAGES" not provided in user_incar_settings, in Angstrom.
    """
    required_params = ["label"]
    optional_params = ["d_img"]

    def run_task(self, fw_spec):

        label = self["label"]
        assert label in ["parent", "ep0", "ep1"] or "neb" in label, "Unknown label!"

        d_img = float(self.get("d_img", 0.7))  # Angstrom
        wf_name = fw_spec["wf_name"]
        src_dir = os.path.abspath(".")
        dest_dir = os.path.join(fw_spec["_fw_env"]["run_dest_root"], wf_name, label)
        shutil.copytree(src_dir, dest_dir)

        # Update fw_spec based on the type of calculations.
        if "neb" in label:
            # Update all relaxed images.
            subs = glob.glob("[0-2][0-9]")
            nimages = len(subs)
            concar_list = ["{:02}/CONTCAR".format(i) for i in range(nimages)[1:-1]]
            images = [Structure.from_file(contcar) for contcar in concar_list]

            # Update the two ending "images".
            images.insert(0, Structure.from_file("00/POSCAR"))
            images.append(Structure.from_file("{:02}/POSCAR".format(nimages - 1)))
            images = [s.as_dict() for s in images]
            neb = fw_spec.get("neb")
            neb.append(images)
            update_spec = {"neb": neb, "_queueadapter": {"nnodes": str(len(images) - 2),
                                                         "nodes": str(len(images) - 2)}}
            # Use neb walltime if it is in fw_spec
            if fw_spec["neb_walltime"] is not None:
                update_spec["_queueadapter"].update({"walltime": fw_spec.get("neb_walltime")})

        elif label in ["ep0", "ep1"]:
            # Update relaxed endpoint structures.
            file = glob.glob("CONTCAR*")[0]
            ep = Structure.from_file(file, False)  # One endpoint

            if fw_spec.get("incar_images"):  # "incar_images": pre-defined image number.
                update_spec = {label: ep.as_dict(),
                               "_queueadapter": {"nnodes": fw_spec["incar_images"],
                                                 "nodes": fw_spec["incar_images"]}}
            else:
                # Calculate number of images if "IMAGES" tag is not provided.
                index = int(label[-1])
                ep_1_dict = fw_spec.get("ep{}".format(1 - index))  # Another endpoint
                try:
                    ep_1 = Structure.from_dict(ep_1_dict)
                except:
                    ep_1 = ep_1_dict

                max_dist = max(get_endpoint_dist(ep, ep_1))
                nimages = round(max_dist / d_img) or 1
                update_spec = {label: ep, "_queueadapter": {"nnodes": int(nimages),
                                                            "nodes": int(nimages)}}
            # Use neb walltime if it is in fw_spec
            if fw_spec["neb_walltime"] is not None:
                update_spec["_queueadapter"].update({"walltime": fw_spec.get("neb_walltime")})

        else:  # label == "parent"
            f = glob.glob("CONTCAR*")[0]
            s = Structure.from_file(f, False)
            ep0, ep1 = get_endpoints_from_index(s, fw_spec["site_indices"])

            update_spec = {"parent": s.as_dict(), "ep0": ep0.as_dict(), "ep1": ep1.as_dict()}

        # Clear current directory.
        for d in os.listdir(src_dir):
            try:
                os.remove(os.path.join(src_dir, d))
            except:
                shutil.rmtree(os.path.join(src_dir, d))

        return FWAction(update_spec=update_spec)


@explicit_serialize
class RunNEBVaspFake(FiretaskBase):
    """
     Vasp Emulator for NEB, which has a different file arrangement. Similar to RunVaspFake class.

     Required params:
         ref_dir (string): Path to reference vasp run directory with input files in the folder named
            'inputs' and output files in the folder named 'outputs'.

     Optional params:
         params_to_check (list): optional list of incar parameters to check.
     """
    required_params = ["ref_dir"]
    optional_params = ["params_to_check"]

    def run_task(self, fw_spec):
        self._get_params()
        self._verify_inputs()
        self._clear_inputs()
        self._generate_outputs()

    def _get_params(self):
        """Define some convenient variables."""

        self.user_dir = os.getcwd()
        self.ref_dir_input = os.path.join(self["ref_dir"], "inputs")
        self.ref_dir_output = os.path.join(self["ref_dir"], "outputs")

        user_sdir = glob.glob(os.path.join(os.getcwd(), "[0-9][0-9]"))
        ref_sdir_input = glob.glob(os.path.join(self["ref_dir"], "inputs", "[0-9][0-9]"))
        ref_sdir_output = glob.glob(os.path.join(self["ref_dir"], "outputs", "[0-9][0-9]"))
        user_sdir.sort()
        ref_sdir_input.sort()
        ref_sdir_output.sort()

        # Check sub-folders consistence.
        if len(user_sdir) != len(ref_sdir_input):
            raise ValueError("Sub-folder numbers are inconsistent! "
                             "Paths are:\n{}\n{}".format(self.user_dir, self.ref_dir_input))
        self.user_sdir = user_sdir
        self.ref_sdir_input = ref_sdir_input
        self.ref_sdir_output = ref_sdir_output

    def _verify_inputs(self):
        """Validation of input files under user NEB directory."""
        user_incar = Incar.from_file(os.path.join(self.user_dir, "INCAR"))
        ref_incar = Incar.from_file(os.path.join(self.ref_dir_input, "INCAR"))

        # Check INCAR
        params_to_check = self.get("params_to_check", [])
        defaults = {"ICHAIN": 0, "LCLIMB": True}
        for p in params_to_check:
            if user_incar.get(p, defaults.get(p)) != ref_incar.get(p, defaults.get(p)):
                raise ValueError("INCAR value of {} is inconsistent!".format(p))

        # Check KPOINTS
        user_kpoints = Kpoints.from_file(os.path.join(self.user_dir, "KPOINTS"))
        ref_kpoints = Kpoints.from_file(os.path.join(self.ref_dir_input, "KPOINTS"))
        if user_kpoints.style != ref_kpoints.style or user_kpoints.num_kpts != ref_kpoints.num_kpts:
            raise ValueError("KPOINT files are inconsistent! "
                             "Paths are:\n{}\n{} with kpts = {} {}".format(
                self.user_dir, self.ref_dir_input, user_kpoints, ref_kpoints))

        # Check POTCAR
        user_potcar = Potcar.from_file(os.path.join(self.user_dir, "POTCAR"))
        ref_potcar = Potcar.from_file(os.path.join(self.ref_dir_input, "POTCAR"))
        if user_potcar.symbols != ref_potcar.symbols:
            raise ValueError("POTCAR files are inconsistent! "
                             "Paths are:\n{}\n{}".format(self.user_dir, self.ref_dir_input))

        # Check POSCARs
        for u, r in zip(self.user_sdir, self.ref_sdir_input):
            user_poscar = Poscar.from_file(os.path.join(u, "POSCAR"))
            ref_poscar = Poscar.from_file(os.path.join(r, "POSCAR"))
            if user_poscar.natoms != ref_poscar.natoms or \
                            user_poscar.site_symbols != ref_poscar.site_symbols:
                raise ValueError("POSCAR files are inconsistent! Paths are:\n{}\n{}".format(u, r))

    def _clear_inputs(self):
        """
        Remove all input files from user NEB directory.
        """
        # Clear neb directory
        for x in VASP_NEB_OUTPUT_FILES:
            p = os.path.join(os.getcwd(), x)
            if os.path.exists(p):
                os.remove(p)

        # Clear neb sub-directory
        for d in self.user_sdir:
            for x in VASP_NEB_OUTPUT_SUB_FILES:
                p = os.path.join(d, x)
                if os.path.exists(p):
                    os.remove(p)

    def _generate_outputs(self):
        """Copy valid outputs from a reference folder to user NEB directory."""
        # Copy NEB files.
        for file_name in os.listdir(self.ref_dir_output):
            full_file_name = os.path.join(self.ref_dir_output, file_name)
            if os.path.isfile(full_file_name):
                shutil.copy(full_file_name, os.getcwd())

        # Copy NEB sub-files.
        for u_dir, r_dir in zip(self.user_sdir, self.ref_sdir_output):
            for file_name in os.listdir(r_dir):
                full_file_name = os.path.join(r_dir, file_name)
                if os.path.isfile(full_file_name):
                    shutil.copy(full_file_name, u_dir)


@explicit_serialize
class WriteNEBFromImages(FiretaskBase):
    """
    Generate CI-NEB input sets using given images and endpoints structures. The structures
    correspond to structures nested in subfolders ("00", "01", "02", etc.)

    Required parameters:
        neb_label (str): "1", "2", etc., labels the running sequence of NEB
    Optional parameters:
        user_incar_settings (dict): Additional INCAR settings.
        user_kpoints_settings (dict): Additional KPOINTS settings.
    """
    required_params = ["neb_label"]
    optional_params = ["user_incar_settings", "user_kpoints_settings"]

    def run_task(self, fw_spec):
        user_incar_settings = self.get("user_incar_settings", {})
        user_kpoints_settings = self.get("user_kpoints_settings", {})
        neb_label = self.get("neb_label")
        assert neb_label.isdigit() and int(neb_label) >= 1
        images = fw_spec["neb"][int(neb_label) - 1]
        try:
            images = [Structure.from_dict(i) for i in images]
        except:
            images = images
        vis = MVLCINEBSet(images, user_incar_settings=user_incar_settings,
                          user_kpoints_settings=user_kpoints_settings)
        vis.write_input(".")


@explicit_serialize
class WriteNEBFromEndpoints(FiretaskBase):
    """
    Generate NEB input sets using endpoint structures, default writes CI-NEB input.
    MVLCINEBSet is the only vasp_input_set supported now.

    The number of images:
        1) search in "user_incar_settings";
        2) otherwise, calculated using "d_img".
    Required parameters:
        user_incar_settings (dict): additional INCAR settings.

    Optional parameters:
        user_kpoints_settings (dict): additional KPOINTS settings.
        sort_tol (float): Distance tolerance (in Angstrom) used to match the atomic indices between
            start and end structures. If it is set 0, then no sorting will be performed.
        d_img (float): distance in Angstrom, used in calculating number of images. Default 0.7
            Angstrom.
        interpolation_type (str): method to do image interpolation from two endpoints.
                            Choose from ["IDPP", "linear"], default "IDPP"
    """
    required_params = ["user_incar_settings"]
    optional_params = ["user_kpoints_settings", "sort_tol", "d_img",
                       "interpolation_type"]

    def run_task(self, fw_spec):
        user_incar_settings = self["user_incar_settings"]
        interpolation_type = self.get("interpolation_type", "IDPP")
        idpp_species = fw_spec.get("idpp_species")
        user_kpoints_settings = self.get("user_kpoints_settings")
        try:
            ep0 = Structure.from_dict(fw_spec["ep0"])
            ep1 = Structure.from_dict(fw_spec["ep1"])
        except:
            ep0 = fw_spec["ep0"]
            ep1 = fw_spec["ep1"]

        # Get number of images.
        nimages = user_incar_settings.get("IMAGES", self._get_nimages(ep0, ep1))
        if interpolation_type == "IDPP":
            from pymatgen_diffusion.neb.pathfinder import IDPPSolver

            obj = IDPPSolver.from_endpoints([ep0, ep1], nimages=nimages)
            images = obj.run(species=idpp_species)
            images_dic_list = [image.as_dict() for image in images]
        elif interpolation_type == "linear":
            images = self._get_images_by_linear_interp(nimages, ep0, ep1)
            images_dic_list = [i.as_dict() for i in images]
        else:
            raise ValueError("The interpolation method must either be 'linear' or 'IDPP'!")

        write = WriteNEBFromImages(neb_label='1', user_incar_settings=user_incar_settings,
                                   user_kpoints_settings=user_kpoints_settings)
        fw_spec["neb"] = [images_dic_list]
        write.run_task(fw_spec=fw_spec)

    def _get_nimages(self, ep0, ep1):
        """
        Calculate the number of images using "d_img", which can be
        overwritten in a optional_params list.

        Returns:
            nimages (int): number of images.
        """
        # Check endpoints consistence.
        assert ep0.atomic_numbers == ep1.atomic_numbers, "Endpoints are inconsistent!"
        d_img = self.get("d_img") or 0.7

        # Assume dilute diffusion
        max_dist = 0
        for s_0, s_1 in zip(ep0, ep1):
            site_dist = ep0.lattice.get_distance_and_image(s_0.frac_coords, s_1.frac_coords)[0]
            if site_dist > max_dist:
                max_dist = site_dist

        # Number of images must more than one.
        nimages = round(max_dist / d_img) or 1

        return nimages

    def _get_images_by_linear_interp(self, nimages, ep0, ep1):
        logger = get_logger(__name__)
        sort_tol = self.get("sort_tol", 0.0)

        try:
            images = ep0.interpolate(ep1, nimages=nimages + 1, autosort_tol=sort_tol)
        except Exception as e:
            if "Unable to reliably match structures " in str(e):
                logger.warn("Auto sorting is turned off because it is unable to match the "
                            "end-point structures!")
                images = ep0.interpolate(ep1, nimages=nimages + 1, autosort_tol=0)
            else:
                raise e

        return images
