# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import glob
import shutil

from pymatgen.core import Structure
from pymatgen.io.vasp import Incar, Kpoints, Poscar, Potcar
from pymatgen.core.surface import SlabGenerator, get_symmetrically_distinct_miller_indices, ReconstructionGenerator
from pymatgen.io.vasp.sets import MVLSlabSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from atomate.vasp.fireworks.core import OptimizeFW
from fireworks.core.firework import FiretaskBase, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize

from atomate.utils.utils import get_logger

"""
Surface workflow firetasks.
"""

__author__ = "Richard Tran"
__email__ = 'rit001@eng.ucsd.edu'

@explicit_serialize
class SurfaceWFGlueTask(FiretaskBase):
    """
    This class generates either a list of child FWs for oriented unit cell
    calculations which subsequently generate child FWs of their own for slab
    calculations or it just directly generates slab calculations as the child
    FWs. In regards to slabs, for now, let's keep it simple and calculate all
    terminations and reconstructions (and for the future, nonstoichiometric,
    tasker 2, etc for multi-element compounds).

    Required params:
        fwname (str): Name of the completed calculation for either a
            conventional unit cell (if this is a node for OUC calculations)
            or OUC (if this is a node for slabs calculations).
        ouc_node (bool): If true, this task generate a list of child slab
            calculation fws from the ouc, otherwise it generates a list of
            child ouc calculation fws from the conventional unit cell

    Optional params:
        surf_qe (QueryEngine): QueryEngine containing all calculations in
            regards to surface properties. If provided, the task will search
            for the final structure in the database using the fwname
        miller_index ([h,k,l]): Miller index of the ouc. Must be provided
            if ouc_node=True
        min_slab_size (in Angstrom):
        min_vac_size (in Angstroms):
    """
    required_params = ["fwname", "ouc_node", "miller_index", "mmi"]
    optional_params = ["db_config", "miller_index", "min_slab_size", "min_vac_size"]

    def run_task(self, fw_spec):

        s = self.get_final_structure()
        min_slab_size = self.get("min_slab_size", 10)
        min_vac_size = self.get("min_vac_size", 10)

        list_of_fws = []
        # generate list of child fw for slabs
        if self.get("ouc_node"):
            # The actual Miller index is not (001), but since the ouc is
            # already oriented in its hkl direction, we want to cleave
            # slabs along that direction which is now the c direction
            slabgen = SlabGenerator(s, (0,0,1), min_slab_size, min_vac_size,
                                    center_slab=True, max_normal_search=max(hkl))
            slabs = slabgen.get_slabs(symmetrize=True)
            for slab in slabs:
                list_of_fws.append(self.slab_optimize_fw(slab))

        # generate list of child fw for oucs
        else:
            symbol = SpacegroupAnalyzer(s).get_space_group_symbol()
            for hkl in get_symmetrically_distinct_miller_indices(s, self.get("mmi")):
                # Check for reconstructions while we're at it. These are
                # constructed from the conventional unit cell. Enumerate
                # through all posisble reconstructions in the archive
                # available for this particular spacegroup
                for rec_name, instructions in reconstructions_archive.items():
                    if "base_reconstruction" in instructions.keys():
                        instructions = reconstructions_archive[instructions["base_reconstruction"]]
                    if instructions["spacegroup"]["symbol"] == symbol:
                        # check if this reconstruction has a max index
                        # equal or less than the given max index
                        if max(instructions["miller_index"]) > max_index:
                            continue
                        recon = ReconstructionGenerator(structure, min_slab_size,
                                                        min_vacuum_size, rec_name)
                        list_of_fws.append(self.slab_optimize_fw(recon.build_slab()))

                ouc = SlabGenerator(s, hkl, 10, 10, max_normal_search=max(hkl)).get_slab().oriented_unit_cell
                list_of_fws.append(self.ouc_optimize_fw(ouc, hkl))

        return FWAction(additions=FWs)

    def get_final_structure(self):

        surf_qe = self.get("surf_qe", None)
        # Queries for the structure in the surface database
        if surf_qe:
            return surf_qe.get_entry({"name": self.get("name")},
                                     inc_structure="Final").structure
        # Load the CONTCAR from a directory with an associated namefile
        else:
            for dir in glob.glob("*"):
                if "FW--{}".format(self.get("name")) in glob.glob("*"):
                    return Structure.from_file("CONTCAR.relax2.gz")

    def slab_optimize_fw(self, slab):

        mvl_slab = MVLSlabSet(slab, k_product=k_product, bulk=False)
        hkl = slab.miller_index
        name = '-%s_slab_k%s_s%ss%s_%s%s%s_shift%s' % (mpid, k_product, min_slab_size, min_vac_size,
                                                       hkl[0], hkl[1], hkl[2], slab.shift)
        return OptimizeFW(slab, name=name, vasp_input_set=mvl_slab,
                          ediffg=-0.02, vasp_cmd=vasp_cmd, parent=parent)

    def ouc_optimize_fw(self, ouc, hkl):

        mvl_ouc = MVLSlabSet(ouc, k_product=k_product, bulk=True)
        name = '-%s_bulk_k%s_%s%s%s' % (mpid, k_product, hkl[0], hkl[1], hkl[2])
        return OptimizeFW(ouc, name=name, vasp_input_set=mvl_ouc,
                          ediffg=-0.02, vasp_cmd=vasp_cmd, parent=parent)
