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

    required_params = ['structure_type', "scratch_dir", "k_product", "vasp_cmd"]
    optional_params = ["slab_gen_params", "mpid", "mmi", "db_file", "miller_index"]

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
            slab_gen_params (dict): Parameters for SlabGenerator or generate_all_slabs
            mpid (str): Materials Project ID of the conventional unit cell.
            mmi (int): Max Miller index.
            miller_index ([h, k, l]): Miller index of plane parallel to
                surface (and oriented unit cell).
            cwd (str): Location of directory to operate the
                workflow and store the final outputs
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """

        # get the parameters for SlabGenerator if we are getting a
        # list of Slab calculation FWs or generate_all_slabs if we
        # are getting a list of oritend unit cell calculation FWs
        if self.get("slab_gen_params", None):
            slab_gen_params = self.get("slab_gen_params", None)
        else:
            slab_gen_params = {"min_slab_size": 10, "min_vacuum_size": 10,
                               "symmetrize": True, "center_slab": True}
            slab_gen_params["max_normal_search"] = self.get("mmi") if \
                self.get("mmi") else max(self.get("miller_index"))
            if self.get("structure_type") == "conventional_unit_cell":
                slab_gen_params["include_reconstructions"] = True

        FWs = []
        if self.get("structure_type") == "conventional_unit_cell":
            # Then we create a set of FWs for oriented_unit_cell
            # calculations and reconstructed slabs

            slab_gen_params["structure"] = Structure.from_file("CONTCAR.relax2.gz")

            all_slabs = generate_all_slabs(max_index=self.get("mmi"),
                                           **slab_gen_params)
            # Get all oriented unit cells up ot a max Miller index (mmi)
            miller_list = []
            for slab in all_slabs:
                if slab.reconstruction:
                    FWs.append(self.get_ouc_fw(slab))
                else:
                    # There are several surface terminations for an oriented unit
                    # cell, we only need to calculate the oriented unit cell once
                    if tuple(slab.miller_index) in miller_list:
                        continue
                    else:
                        # build a oriented unit cell fw
                        FWs.append(self.get_ouc_fw(slab))
                        miller_list.append(tuple(slab.miller_index))

        elif self.get("structure_type") == "oriented_unit_cell":

            folder = os.path.basename(os.getcwd())
            # If this is a reconstruction, we need to use the ReconstructionGenerator
            if "_rec_" in folder:
                ouc = Structure.from_file("CONTCAR.relax2.gz")
                # ReconstructionGenerator only works on the conventional ucell
                ucell = SpacegroupAnalyzer(ouc).get_conventional_standard_structure()
                # Get the name of the reocnstruction
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
                FWs.append(self.get_slab_fw(rec.build_slab(), slab_gen_params))

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
                    slab.miller_index = self.get("miller_index")
                    FWs.append(self.get_slab_fw(slab, slab_gen_params))

        return FWAction(additions=FWs)

    def get_ouc_fw(self, slab):
        """
        Return a SurfCalcOptimizer FW for an oriented unit cell.

        Args:
            slab (Slab): Slab object containing various
                attributes related to the slab
        """

        from atomate.vasp.fireworks.core import SurfCalcOptimizer

        return SurfCalcOptimizer(slab.oriented_unit_cell, self.get("scratch_dir"),
                                 self.get("k_product"), self.get("db_file"),
                                 self.get("vasp_cmd"), "oriented_unit_cell",
                                 cwd=os.getcwd(),
                                 reconstruction=slab.reconstruction,
                                 miller_index=slab.miller_index,
                                 scale_factor=slab.scale_factor,
                                 mpid=self.get("mpid", "--"))

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

        return SurfCalcOptimizer(slab, self.get("scratch_dir"),
                                 self.get("k_product"), self.get("db_file"),
                                 self.get("vasp_cmd"), "slab_cell",
                                 cwd=os.getcwd(),
                                 miller_index=slab.miller_index,
                                 scale_factor=slab.scale_factor,
                                 ouc=slab.oriented_unit_cell, shift=slab.shift,
                                 ssize=slab_gen_params["min_slab_size"],
                                 vsize=slab_gen_params["min_vacuum_size"],
                                 reconstruction=slab.reconstruction,
                                 mpid=self.get("mpid", "--"))