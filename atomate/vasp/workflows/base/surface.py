# coding: utf-8


## for Surface Energy Calculation
from __future__ import division, unicode_literals

__author__ = "Richard Tran"
__version__ = "0.2"
__email__ = "rit001@eng.ucsd.edu"
__date__ = "11/30/17"

import json
import os

from atomate.vasp.firetasks.run_calc import RunVaspCustodian
from atomate.common.firetasks.glue_tasks import CreateFolder, RenameFile
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSet
from atomate.vasp.firetasks.parse_outputs import VaspToDb

from pymatgen.core.surface import generate_all_slabs, Structure, SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.sets import MVLSlabSet

from fireworks.core.firework import Firework, Workflow, FiretaskBase, FWAction
from fireworks import explicit_serialize


class SurfaceWorkflowManager(object):
    """
        Initializes the workflow manager by taking in a list of
        formulas/mpids or a dictionary with the formula as the key referring
        to a list of miller indices.
    """

    def __init__(self, db_file, cwd=os.getcwd(),
                 scratch_dir="", k_product=50, vasp_cmd="vasp"):

        with open(db_file) as data_file:
            dbconfig = json.load(data_file)

        self.k_product = k_product
        self.vasp_cmd = vasp_cmd
        self.dbconfig = dbconfig
        self.db_file = db_file
        self.scratch_dir = scratch_dir
        self.cwd = cwd

    def from_conventional_unit_cell(self, structure, mmi, mpid="--"):

        return Workflow([SurfCalcOptimizer(structure, self.scratch_dir,
                                           self.k_product, self.db_file,
                                           self.vasp_cmd, "conventional_unit_cell",
                                           mmi=mmi, cwd=self.cwd, mpid=mpid)])

    def from_oriented_unit_cell(self, structure, miller_index, scale_factor, mpid="--"):

        return Workflow([SurfCalcOptimizer(structure, self.scratch_dir, self.k_product,
                                           self.db_file, self.vasp_cmd, "oriented_unit_cell",
                                           miller_index=miller_index, scale_factor=scale_factor,
                                           cwd=self.cwd, mpid=mpid, **kwargs)])

    def from_slab_cell(self, structure, miller_index, shift, scale_factor,
                       ouc, ssize, vsize, reconstruction=None, mpid="--"):

        return Workflow([SurfCalcOptimizer(structure, self.scratch_dir, self.k_product,
                                           self.db_file, self.vasp_cmd, "slab_cell",
                                           miller_index=miller_index, ouc=ouc, shift=shift,
                                           scale_factor=scale_factor, ssize=ssize, vsize=vsize,
                                           reconstruction=reconstruction, cwd=self.cwd,
                                           mpid=mpid, **kwargs)])


class SurfCalcOptimizer(Firework):
    def __init__(self, structure, scratch_dir, k_product, db_file,
                 vasp_cmd, structure_type, miller_index=[], scale_factor=[],
                 vsize=None, mmi=None, ouc=None, shift=None, ssize=None,
                 reconstruction=None, cwd=os.getcwd(), mpid="--", **kwargs):
        """
        Customized FW similar to OptimizeFW.

        Args:
            scratch_dir: (str) - if specified, uses this directory as the root
                scratch dir. Supports env_chk.
            k_product (int): Default to 50, kpoint number * length for a & b
                directions, also for c direction in bulk calculations.
            db_file (str): FULL path to file containing the database credentials.
                Supports env_chk.
            vasp_cmd (str): Command to run vasp.
            structure (Structure): Input structure. Need to input either the
                structure or an mpid
            mpid (str): Unique Materials Project structure ID of the unit cell.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """

        self.structure = structure
        self.el = self.structure[0].species_string
        self.structure_type = structure_type
        self.k_product = k_product
        self.mpid = mpid
        self.hkl = miller_index
        self.reconstruction = reconstruction
        self.ssize = ssize
        self.vsize = vsize
        self.shift = shift
        self.sg = SpacegroupAnalyzer(ouc if structure_type ==
                                            "slab_cell" else structure)
        self.scale_factor = scale_factor
        self.ouc = ouc
        self.cwd = cwd
        self.vasp_cmd = vasp_cmd
        self.scratch_dir = scratch_dir
        self.db_file = db_file
        self.mmi = mmi

        super(SurfCalcOptimizer, self).__init__(self.get_tasks, name=self.get_name, **kwargs)

    @property
    def get_input_set(self):

        if self.structure_type != "slab_cell":
            return MVLSlabSet(self.structure, bulk=True,
                              k_product=self.k_product)
        else:
            return MVLSlabSet(self.structure, bulk=True,
                              k_product=self.k_product)

    @property
    def get_name(self):
        if self.structure_type == "conventional_unit_cell":
            return "%s_%s_conventional_unit_cell_k%s" % \
                   (self.el, self.mpid, self.k_product)
        elif self.structure_type == "oriented_unit_cell":
            return "%s_%s_bulk_k%s_%s%s%s" % \
                   (self.el, self.mpid, self.k_product,
                    self.hkl[0], self.hkl[1], self.hkl[2])
        elif self.structure_type == "slab_cell":
            if not self.reconstruction:
                return "%s_%s_slab_k%s_s%sv%s_%s%s%s_shift%s" % \
                       (self.el, self.mpid, self.k_product, self.ssize, self.vsize,
                        self.hkl[0], self.hkl[1], self.hkl[2], self.shift)
            else:
                return "%s_%s_slab_k%s_s%sv%s_%s" % (self.el, self.mpid, \
                                                     self.k_product, self.ssize,
                                                     self.vsize, self.reconstruction)

    @property
    def get_additional_fields(self):

        additional_fields = {"structure_type": self.structure_type,
                             "calculation_name": self.get_name,
                             "conventional_spacegroup": \
                                 {"symbol": self.sg.get_space_group_symbol(),
                                  "number": self.sg.get_space_group_number()},
                             "initial_structure": self.structure.as_dict(),
                             "material_id": self.mpid}

        if self.structure_type != "conventional_unit_cell":
            additional_fields.update({"miller_index": self.hkl,
                                      "scale_factor": self.scale_factor})

        if self.structure_type == "slab_cell":
            additional_fields.update({"oriented_unit_cell": self.ouc.as_dict(),
                                      "slab_size":self.ssize, "shift": self.shift,
                                      "vac_size": self.vsize,
                                      "reconstruction": self.reconstruction})

        return additional_fields

    @property
    def get_tasks(self):

        tasks = [CreateFolder(folder_name=os.path.join(self.cwd, self.get_name),
                              change_dir=True, relative_path=True),
                 WriteVaspFromIOSet(structure=self.structure,
                                    vasp_input_set=self.get_input_set),
                 RunVaspCustodian(vasp_cmd=self.vasp_cmd,
                                  scratch_dir=self.scratch_dir,
                                  auto_npar=">>auto_npar<<",
                                  job_type="double_relaxation_run")]

        if self.structure_type == "slab_cell":
            tasks.append(RenameFile(file="LOCPOT.gz", new_name="LOCPOT.relax2.gz"))

        tasks.append(VaspToDb(additional_fields=self.get_additional_fields,
                              db_file=self.db_file))

        if self.structure_type != "slab_cell":
            tasks.append(FacetFWsGeneratorTask(structure_type=self.structure_type,
                                               vasp_cmd=self.vasp_cmd, cwd=self.cwd,
                                               db_file=self.db_file, miller_index=self.hkl,
                                               scratch_dir=self.scratch_dir, mmi=self.mmi,
                                               mpid=self.mpid, k_product=self.k_product))

        return tasks


@explicit_serialize
class FacetFWsGeneratorTask(FiretaskBase):
    """
    Task for generating FWs for oriented unit cell calculations or slab
    calculations. If the initial structure is an oriented unit cell,
    generate slab fws, if its a conventional unit cell, generate oriented
    unit cell and reconstruction calculations.
    """

    required_params = ['structure_type', "scratch_dir", "k_product",
                       "db_file", "vasp_cmd"]
    optional_params = ["slab_gen_params", "mpid", "mmi", "miller_index", "cwd"]

    def run_task(self, fw_spec):

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
            miller_list = []
            for slab in all_slabs:
                if slab.reconstruction:
                    FWs.append(self.get_slab_fw(slab, slab_gen_params))
                else:
                    if tuple(slab.miller_index) in miller_list:
                        continue
                    else:
                        # build a oriented unit cell fw
                        FWs.append(self.get_ouc_fw(slab))

        elif self.get("structure_type") == "oriented_unit_cell":
            slab_gen_params["initial_structure"] = \
                Structure.from_file("CONTCAR.relax2.gz").as_dict()
            slab_gen_params["miller_index"] = [0,0,1]
            symmetrize = slab_gen_params["symmetrize"]
            del slab_gen_params["symmetrize"]
            slabgen = SlabGenerator(**slab_gen_params)
            for slab in slabgen.get_slabs(symmetrize=symmetrize):
                slab.miller_index = self.get("miller_index")
                FWs.append(self.get_slab_fw(slab, slab_gen_params))

        return FWAction(additions=FWs)

    def get_ouc_fw(self, slab):

        return SurfCalcOptimizer(slab.oriented_unit_cell, self.get("scratch_dir"),
                                 self.get("k_product"), self.get("db_file"),
                                 self.get("vasp_cmd"), "oriented_unit_cell",
                                 miller_index=slab.miller_index,
                                 scale_factor=slab.scale_factor,
                                 mpid=self.get("mpid", "--"))

    def get_slab_fw(self, slab, slab_gen_params):

        return SurfCalcOptimizer(slab.oriented_unit_cell, self.get("scratch_dir"),
                                 self.get("k_product"), self.get("db_file"),
                                 self.get("vasp_cmd"), "slab_cell",
                                 miller_index=slab.miller_index,
                                 scale_factor=slab.scale_factor,
                                 ouc=slab.oriented_unit_cell, shift=slab.shift,
                                 ssize=slab_gen_params["min_slab_size"],
                                 vsize=slab_gen_params["min_vacuum_size"],
                                 reconstruction=slab.reconstruction,
                                 mpid=self.get("mpid", "--"))