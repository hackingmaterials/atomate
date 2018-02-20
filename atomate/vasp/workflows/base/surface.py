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
                 production_mode=False,
                 scratch_dir="", k_product=50, vasp_cmd="vasp"):

        with open(db_file) as data_file:
            dbconfig = json.load(data_file)

        self.k_product = k_product
        self.vasp_cmd = vasp_cmd
        self.dbconfig = dbconfig
        self.db_file = db_file
        self.scratch_dir = scratch_dir
        self.cwd = cwd

    def from_conventional_unit_cell(self, structure, mmi, mpid=None):

        return Workflow([ConvUcellFW(structure, mmi,
                                     self.production_mode, self.scratch_dir,
                                     self.k_product, self.db_file, self.vasp_cmd,
                                     cwd=self.cwd, mpid=mpid)])

    def from_oriented_unit_cell(self, structure, miller_index, mpid=None):

        return Workflow([OUCFW(ouc, miller_index, scale_factor, scratch_dir,
                               k_product, db_file, vasp_cmd, cwd=os.getcwd(),
                               mpid="--", **kwargs)])

    def from_slab_cell(self, structure, miller_index, mpid=None):

        return Workflow([SlabFW(slab, miller_index, scale_factor, ouc, shift,
                                ssize, vsize, scratch_dir, k_product, db_file,
                                vasp_cmd, reconstruction=None, cwd=os.getcwd(),
                                mpid="--", **kwargs)])


class ConvUcellFW(Firework):
    def __init__(self, ucell, mmi, scratch_dir, k_product,
                 db_file, vasp_cmd, cwd=os.getcwd(), mpid="--", **kwargs):
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

        name = "%s_%s_conventional_unit_cell_k%s" % \
               (ucell[0].species_string, mpid, k_product)
        mvl = MVLSlabSet(ucell, k_product=k_product, bulk=True)

        # This will give us our four basic tasks: WriteVaspFromIOSet,
        # RunVaspCustodian, PassCalcLocs and VaspToDB. We can then
        # modify or remove tasks to suit our needs.
        tasks = []
        tasks.append(CreateFolder(folder_name=os.path.join(cwd, name),
                                  change_dir=True, relative_path=True))
        tasks.append(WriteVaspFromIOSet(structure=ucell, vasp_input_set=mvl))
        tasks.append(RunVaspCustodian(vasp_cmd=vasp_cmd,
                                      auto_npar=">>auto_npar<<",
                                      job_type="double_relaxation_run",
                                      scratch_dir=scratch_dir))

        # Add additional fields to distinguish this calculation
        additional_fields = {"structure_type": "conventional_unit_cell",
                             "calculation_name": name}
        additional_fields["conventional_spacegroup"] = \
            SpacegroupAnalyzer(ucell).get_space_group_symbol()
        additional_fields["initial_structure"] = ucell.as_dict()
        additional_fields["material_id"] = mpid
        tasks.append(VaspToDb(additional_fields=additional_fields, db_file=db_file))
        tasks.append(FacetFWsGeneratorTask(structure_type="conventional_unit_cell", mmi=mmi,
                                           vasp_cmd=vasp_cmd, mpid=mpid, db_file=db_file,
                                           scratch_dir=scratch_dir, k_product=k_product, cwd=cwd))

        super(ConvUcellFW, self).__init__(tasks, name=name, **kwargs)

class OUCFW(Firework):
    def __init__(self, ouc, miller_index,
                 scale_factor, scratch_dir, k_product,
                 db_file, vasp_cmd, cwd=os.getcwd(), mpid="--", **kwargs):
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

        name = "%s_%s_bulk_k%s_%s%s%s" % \
               (ouc[0].species_string, mpid, k_product,
                miller_index[0], miller_index[1], miller_index[2])
        mvl = MVLSlabSet(ouc, k_product=k_product, bulk=True)

        # This will give us our four basic tasks: WriteVaspFromIOSet,
        # RunVaspCustodian, PassCalcLocs and VaspToDB. We can then
        # modify or remove tasks to suit our needs.
        tasks = []
        tasks.append(CreateFolder(folder_name=os.path.join(cwd, name),
                                  change_dir=True, relative_path=True))
        tasks.append(WriteVaspFromIOSet(structure=ouc, vasp_input_set=mvl))
        tasks.append(RunVaspCustodian(vasp_cmd=vasp_cmd,
                                      auto_npar=">>auto_npar<<",
                                      job_type="double_relaxation_run",
                                      scratch_dir=scratch_dir))

        # Add additional fields to distinguish this calculation
        additional_fields = {"structure_type": "oriented_unit_cell",
                             "calculation_name": name}
        additional_fields["conventional_spacegroup"] = \
            SpacegroupAnalyzer(ouc).get_space_group_symbol()
        additional_fields["initial_structure"] = ouc.as_dict()
        additional_fields["material_id"] = mpid
        additional_fields["miller_index"] = miller_index
        additional_fields["scale_factor"] = scale_factor
        tasks.append(VaspToDb(additional_fields=additional_fields, db_file=db_file))
        tasks.append(FacetFWsGeneratorTask(structure_type="oriented_unit_cell",
                                           vasp_cmd=vasp_cmd, mpid=mpid, db_file=db_file,
                                           scratch_dir=scratch_dir, k_product=k_product,
                                           miller_index=miller_index, cwd=cwd))

        super(OUCFW, self).__init__(tasks, name=name, **kwargs)

class SlabFW(Firework):
    def __init__(self, slab,
                 miller_index, scale_factor, ouc, shift, ssize, vsize,
                 scratch_dir, k_product, db_file,
                 vasp_cmd, reconstruction=None, cwd=os.getcwd(),
                 mpid="--", **kwargs):
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

        if not reconstruction:
            name = "%s_%s_slab_k%s_s%sv%s_%s%s%s_shift%s" % \
                   (slab[0].species_string, mpid, k_product, ssize, vsize,
                    miller_index[0], miller_index[1], miller_index[2], shift)
        else:
            name = "%s_%s_slab_k%s_s%sv%s_%s" % (slab[0].species_string, mpid, \
                                                 k_product, ssize, vsize, reconstruction)

        mvl = MVLSlabSet(slab, k_product=k_product, bulk=False, get_wf=True)

        # This will give us our four basic tasks: WriteVaspFromIOSet,
        # RunVaspCustodian, PassCalcLocs and VaspToDB. We can then
        # modify or remove tasks to suit our needs.
        tasks = []
        tasks.append(CreateFolder(folder_name=os.path.join(cwd, name),
                                  change_dir=True, relative_path=True))
        tasks.append(WriteVaspFromIOSet(structure=slab, vasp_input_set=mvl))
        tasks.append(RunVaspCustodian(vasp_cmd=vasp_cmd,
                                      auto_npar=">>auto_npar<<",
                                      job_type="double_relaxation_run",
                                      scratch_dir=scratch_dir))

        # Add additional fields to distinguish this calculation
        additional_fields = {"structure_type": "slab_cell",
                             "calculation_name": name}
        sg = SpacegroupAnalyzer(ouc)
        additional_fields["conventional_spacegroup"] = \
            {"symbol": sg.get_space_group_symbol(), "number": sg.get_space_group_number()}
        additional_fields["initial_structure"] = slab.as_dict()
        additional_fields["material_id"] = mpid
        additional_fields["miller_index"] = miller_index
        additional_fields["scale_factor"] = scale_factor
        additional_fields["oriented_unit_cell"] = ouc.as_dict()
        additional_fields["slab_size"] = ssize
        additional_fields["vac_size"] = vsize
        additional_fields["shift"] = shift
        additional_fields["reconstruction"] = reconstruction
        tasks.append(RenameFile(file="LOCPOT.gz", new_name="LOCPOT.relax2.gz"))
        tasks.append(VaspToDb(additional_fields=additional_fields, db_file=db_file))
        super(SlabFW, self).__init__(tasks, name=name, **kwargs)


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

        return OUCFW(slab.oriented_unit_cell, slab.miller_index,
                     slab.scale_factor,
                     self.get("scratch_dir"), self.get("k_product"),
                     self.get("db_file"), self.get("vasp_cmd"),
                     mpid=self.get("mpid", "--"))

    def get_slab_fw(self, slab, slab_gen_params):

        return SlabFW(slab, slab_gen_params,
                      slab.miller_index, slab.scale_factor,
                      slab.oriented_unit_cell, slab.shift,
                      slab_gen_params["min_slab_size"],
                      slab_gen_params["min_vacuum_size"],
                      self.get("scratch_dir"), self.get("k_product"),
                      self.get("db_file"), self.get("vasp_cmd"),
                      reconstruction=slab.reconstruction,
                      cwd=self.get("cwd", os.getcwd()),
                      mpid=self.get("mpid", "--"))

