# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os

from pymatgen.core import Structure
from pymatgen.core.surface import SlabGenerator, ReconstructionGenerator, generate_all_slabs
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from fireworks.core.firework import FiretaskBase, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize


"""
Surface workflow firetasks.
"""

__author__ = "Richard Tran"
__email__ = 'rit001@eng.ucsd.edu'


@explicit_serialize
class FacetFWsGeneratorTask(FiretaskBase):
    """
    Task for generating FWs for oriented unit cell calculations or slab
    calculations. If the initial structure is an oriented unit cell,
    generate slab fws, if its a conventional unit cell, generate oriented
    unit cell and reconstruction calculations.
    """

    required_params = ['structure_type', "k_product", "vasp_cmd", "scratch_dir"]
    optional_params = ["slab_gen_params", "naming_tag", "max_index",
                       "db_file", "miller_index", "run_dir"]

    def run_task(self, fw_spec):
        """
        Generates a list of oriented unit cell FWs or
        slab cell FWs depending on the structure type.

        Args:
            structure_type (str): Type of structure. Options are:
                conventional_unit_cell, oriented_unit_cell and slab_cell
            scratch_dir (str): - if specified, uses this directory as the root
                scratch dir. Supports env_chk.
            k_product (int): Default to 50, kpoint number * length for a & b
                directions, also for c direction in bulk calculations.
            db_file (str): FULL path to file containing the database credentials.
                Supports env_chk.
            vasp_cmd (str): Command to run vasp.
            slab_gen_params (dict): Parameters for SlabGenerator or
                generate_all_slabs. Defaults to a minimum slab/vacuum
                size of 10A and with a max_normal_search equal to the max Miller index
            naming_tag (str): Naming tag associated with the calculation. Defaults to "--".
            max_index (int): Max Miller index. This is needed if you are starting
                from a conventional_unit_cell.
            miller_index ([h, k, l]): Miller index of plane parallel to
                surface (and oriented unit cell). This is needed if you
                are starting from an oriented_unit_cell or slab_cell.
            run_dir (str): Location of directory to operate the
                workflow and store the final outputs. Defaults
                to your current working directory.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """

        structure_type = self["structure_type"]

        # get the parameters for SlabGenerator if we are getting a
        # list of Slab calculation FWs or generate_all_slabs if we
        # are getting a list of oritend unit cell calculation FWs
        if self.get("slab_gen_params"):
            slab_gen_params = self["slab_gen_params"]
        else:
            slab_gen_params = {"min_slab_size": 10, "min_vacuum_size": 10,
                               "symmetrize": True, "center_slab": True}
            slab_gen_params["max_normal_search"] = self["max_index"] if \
                self.get("max_index") else max(self["miller_index"])
            if structure_type == "conventional_unit_cell":
                slab_gen_params["include_reconstructions"] = True

        FWs = []
        if structure_type == "conventional_unit_cell":
            # Then we create a set of FWs for oriented_unit_cell
            # calculations and reconstructed slabs
            slab_gen_params["structure"] = Structure.from_file("CONTCAR.relax2.gz")

            all_slabs = generate_all_slabs(max_index=self["max_index"],
                                           **slab_gen_params)
            # Get all oriented unit cells up to a max Miller index (max_index)
            miller_list, recon_unit_vects = [], []
            for slab in all_slabs:
                if slab.reconstruction:
                    # Some reconstructions may be based on the same type on oriented
                    # unit cell, avoid duplicate calculations of oriented unit cell
                    m = tuple([tuple(v) for v in slab.recon_trans_matrix])
                    if m in recon_unit_vects:
                        continue
                    else:
                        recon_unit_vects.append(m)
                        FWs.append(self.get_oriented_ucell_fw(slab))
                else:
                    # There are several surface terminations for an oriented unit
                    # cell, we only need to calculate the oriented unit cell once
                    if tuple(slab.miller_index) in miller_list:
                        continue
                    else:
                        # build a oriented unit cell fw
                        FWs.append(self.get_oriented_ucell_fw(slab))
                        miller_list.append(tuple(slab.miller_index))

        elif structure_type == "oriented_unit_cell":
            folder = os.path.basename(os.getcwd())
            # If this is a reconstruction, we need to use the ReconstructionGenerator
            if "_rec_" in folder:
                oriented_ucell = Structure.from_file("CONTCAR.relax2.gz")
                # ReconstructionGenerator only works on the conventional ucell
                ucell = SpacegroupAnalyzer(oriented_ucell).\
                    get_conventional_standard_structure()
                # Get the name of the reconstruction
                ns, n = "", 0
                for s in folder:
                    ns+=s
                    if s == "_":
                        n+=1
                    if n == 5:
                        break
                rec = ReconstructionGenerator(ucell, slab_gen_params["min_slab_size"],
                                              slab_gen_params["min_vacuum_size"],
                                              reconstruction_name=folder[len(ns):])
                for slab in rec.build_slabs():
                    FWs.append(self.get_slab_fw(slab, slab_gen_params))

            else:
                # Get the list of FWs for the various terminations
                # for a slab based on the oriented unit cell
                slab_gen_params["initial_structure"] = \
                    Structure.from_file("CONTCAR.relax2.gz")
                slab_gen_params["miller_index"] = [0,0,1]
                symmetrize = slab_gen_params["symmetrize"]
                del slab_gen_params["symmetrize"]
                slabgen = SlabGenerator(**slab_gen_params)
                for slab in slabgen.get_slabs(symmetrize=symmetrize):
                    slab.miller_index = self["miller_index"]
                    FWs.append(self.get_slab_fw(slab, slab_gen_params))

        return FWAction(additions=FWs)

    def get_oriented_ucell_fw(self, slab):
        """
        Return a SurfCalcOptimizer FW for an oriented unit cell.

        Args:
            slab (Slab): Slab object containing various
                attributes related to the slab
        """

        from atomate.vasp.fireworks.core import SurfCalcOptimizer
        return SurfCalcOptimizer(slab.oriented_unit_cell, self["scratch_dir"],
                                 self["k_product"], self["vasp_cmd"],
                                 "oriented_unit_cell", self.get("run_dir", os.getcwd()),
                                 reconstruction=slab.reconstruction,
                                 miller_index=slab.miller_index,
                                 db_file=self.get("db_file"),
                                 scale_factor=slab.scale_factor,
                                 naming_tag=self.get("naming_tag", "--"))

    def get_slab_fw(self, slab, slab_gen_params):
        """
        Return a SurfCalcOptimizer FW for a Slab cell.

        Args:
            slab (Slab): Slab object containing various
                attributes related to the slab
            slab_gen_params (dict): Parameters for SlabGenerator
                or generate_all_slabs
        """

        from atomate.vasp.fireworks.core import SurfCalcOptimizer

        return SurfCalcOptimizer(slab, self["scratch_dir"],
                                 self["k_product"], self["vasp_cmd"],
                                 "slab_cell", self.get("run_dir", os.getcwd()),
                                 miller_index=slab.miller_index,
                                 db_file=self.get("db_file"),
                                 scale_factor=slab.scale_factor,
                                 oriented_ucell=slab.oriented_unit_cell,
                                 shift=slab.shift,
                                 min_slab_size=slab_gen_params["min_slab_size"],
                                 min_vac_size=slab_gen_params["min_vacuum_size"],
                                 reconstruction=slab.reconstruction,
                                 naming_tag=self.get("naming_tag", "--"))