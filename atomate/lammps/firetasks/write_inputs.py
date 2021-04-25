# coding: utf-8


"""
This module defines firetasks for writing LAMMPS input files (data file and the control
parameters file)
"""

from pymatgen.io.lammps.inputs import LammpsRun, LammpsData
from pymatgen.core import Molecule

from fireworks import FiretaskBase, explicit_serialize

__author__ = 'Kiran Mathew, Brandon Wood, Eric Sivonxay'
__email__ = "kmathew@lbl.gov, b.wood@berkeley.edu, esivonxay@lbl.gov"

@explicit_serialize
class WriteInputFromTemplate(FiretaskBase):
    """

    """

    required_params = ["settings", "lammps_data", "input_filename"]

    optional_params = ["script_template", "template_dir"]

    def run_task(self, fw_spec):
        script_template = self.get("script_template", None)
        settings = self["settings"]
        lammps_data = self["lammps_data"]
        if script_template is None:
            template_dir = self.get("template_dir", None)
            if template_dir is None:
                raise ValueError("no template supplied, defaulting to default md script")
            else:
                with open(template_dir, 'r') as f:
                    script_template = f.read()


        lr = LammpsRun(script_template, settings, lammps_data, self["input_filename"])
        lr.write_inputs('./')


@explicit_serialize
class WriteInputFromIOSet(FiretaskBase):
    """
    Writes LAMMPS Input files(data file and the control parameters file) from DictLammpsInput.

    required_params:
        lammps_input_set (LammpsInputSet)
        input_file (string): name of the file to which the input params will be written

    optional_params:
        data_filename (string): if specified the data file will be renamed
    """

    required_params = ["lammps_input_set", "input_filename"]

    optional_params = ["data_filename"]

    def run_task(self, fw_spec):

        lammps_input_set = self["lammps_input_set"]
        input_filename = self["input_filename"]
        data_filename = self.get("data_filename", None)
        if isinstance(lammps_input_set.lammps_data, dict):
            lammps_input_set.lammps_data = LammpsData.from_dict(lammps_input_set.lammps_data)
        lammps_input_set.write_input(input_filename, data_filename)


@explicit_serialize
class WriteInputFromForceFieldAndTopology(FiretaskBase):
    # Note: This function only works versions of pymatgen from ~2017 and prior
    # TODO: Fix this to work with the current Lammps IO in pymatgen.

    required_params = ["input_file", "final_molecule", "constituent_molecules", "mols_number",
                       "box_size", "forcefield", "topologies", "input_filename"]

    optional_params = ["data_filename", "user_settings", "ff_site_property"]

    def run_task(self, fw_spec):

        molecules = self["constituent_molecules"]
        mols_number = self["mols_number"]
        input_filename = self["input_filename"]
        forcefield = self["forcefield"]
        topologies = self["topologies"]

        user_settings = self.get("user_settings", {})
        data_filename = self.get("data_filename", user_settings.get("data_file", "lammps.data"))
        final_molecule = self["final_molecule"]

        # if the final molecule was generated using packmol
        if fw_spec.get("packed_mol", None):
            final_molecule = fw_spec["packed_mol"]
        elif isinstance(final_molecule, str):
            final_molecule = Molecule.from_file(final_molecule)

        #molecules, mols_number, final_molecule
        lammps_ff_data = LammpsData.from_ff_and_topologies(forcefield, topologies, self["box_size"])

        lammps_input_set = LammpsInputSet.from_file("ff-inputset", self["input_file"],
                                                    user_settings=user_settings,
                                                    lammps_data=lammps_ff_data,
                                                    data_filename=data_filename,
                                                    is_forcefield=True)

        lammps_input_set.write_input(input_filename, data_filename)
