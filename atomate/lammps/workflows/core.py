# coding: utf-8


"""
This module defines functions that yield lammps workflows
"""

from fireworks import Workflow

# from pymatgen.io.lammps.sets import LammpsInputSet
# from pymatgen.io.lammps.data import Topology

from atomate.lammps.fireworks.core import LammpsFW, PackmolFW, LammpsForceFieldFW

__author__ = 'Kiran Mathew'
__email__ = "kmathew@lbl.gov"


def get_wf_basic(input_file, user_settings, lammps_data=None, input_filename="lammps.in",
                 is_forcefield=False, lammps_cmd="lmp_serial", dump_filenames=None, db_file=None,
                 name="LAMMPS Wflow"):
    """
    Returns basic lammps workflow. This is more useful if the input_file is template file with
    the corresponding settings defined in user_settings

    Args:
        input_file (str): path to lammps input file.
            Note: It could be a template file too, then the user_settings must be set.
        user_settings ([dict] or dict): list of settings dict. if the input_file is a tempalte file
            then each dict contains the key value pairs for the template file.
        lammps_data (string/LammpsData/LammpsForceFieldData): path to the data file or
            an appropriate object.
        input_filename (string): input file name. This is the name of the input file passed to the
            lammps binary.
        is_forcefield (bool): whether the data file has forcefield and topology info in it.
            This is required only if lammps_data is a path to the data file instead of a data object.
        lammps_cmd (string): lammps command to run (skip the input file).
        dump_filenames ([str]): list of dump file names
        db_file (string): path to the db file.
        name (str): workflow name

    Returns:
        Workflow
    """
    wf_name = name
    user_settings = user_settings or {}
    user_settings = user_settings if isinstance(user_settings, list) else [user_settings]

    fws = []
    for settings in user_settings:

        data_filename = settings.get("data_file", "lammps.data")

        if "log_file" not in settings:
            settings["log_file"] = "log.lammps"
        log_filename = settings["log_file"]

        lammps_input_set = LammpsInputSet.from_file(wf_name, input_file,
                                                    user_settings=settings, lammps_data=lammps_data,
                                                    data_filename=data_filename)

        fws.append(
            LammpsFW(lammps_input_set=lammps_input_set, input_filename=input_filename,
                     data_filename=data_filename, lammps_cmd=lammps_cmd, db_file=db_file,
                     log_filename=log_filename, dump_filename=dump_filenames)
        )

    return Workflow(fws, name=name)


def get_packmol_wf(input_file, user_settings, constituent_molecules, packing_config, forcefield,
                   final_box_size, topologies=None, ff_site_property=None, tolerance=2.0, filetype="xyz",
                   control_params=None, lammps_cmd="lmp_serial", packmol_cmd="packmol",
                   dump_filenames=None, db_file=None, name="Packmol Lammps Wflow"):
    """
    Returns workflow that uses Packmol to pack the constituent molecules into the given
    configuration and then run lammps on the final packed molecule for the given list of
    user_settings.

    Args:
        input_file (str):  path to lammps input(or template) file.
        user_settings ([dict] or dict): list of settings dict. if the input_file is a tempalte file
            then each dict contains the key value pairs for the template file.
        constituent_molecules ([Molecules]): list of pymatgen Molecule objects
        packing_config ([dict]): list of configuration dictionaries, one for each constituent molecule.
        forcefield (ForceField): pymatgen.io.lammps.forcefield.ForceField object
        final_box_size ([list]): list of list of low and high values for each dimension [[xlow, xhigh], ...]
        topologies ([Topology]): list of Topology objects. If not given, will be set from the
            topology of the constituent molecules.
        ff_site_property (str): the name of the site property used for forcefield mapping
        tolerance (float): packmol tolerance
        filetype (str): packmol i/o file type.
        control_params (dict): packmol control params
        lammps_cmd (string): lammps command to run (skip the input file).
        packmol_cmd (string): path to packmol bin
        dump_filenames ([str]): list of dump file names
        db_file (string): path to the db file.
        name (str): workflow name

    Returns:
        Workflow
    """

    user_settings = user_settings if isinstance(user_settings, list) else [user_settings]

    packmol_output_file = "packed_mol.{}".format(filetype)
    mols_number = [mol_config["number"] for mol_config in packing_config]

    topologies = topologies or []
    # if not given then get the topology from the constituent molecules.
    if not topologies:
        topologies = [Topology.from_molecule(mol, ff_map=ff_site_property or "ff_map")
                      for mol in constituent_molecules]

    fws = []

    fw_packmol = PackmolFW(constituent_molecules, packing_config, tolerance=tolerance,
                           filetype=filetype, control_params=control_params,
                           copy_to_current_on_exit=True, output_file=packmol_output_file,
                           site_property=ff_site_property, packmol_cmd=packmol_cmd)

    fws.append(fw_packmol)

    for setting in user_settings:

        data_filename = setting.get("data_file", "data.lammps")

        if "log_file" not in setting:
            setting["log_file"] = "log.lammps"
        log_filename = setting["log_file"]

        fw_lammps = LammpsForceFieldFW(input_file, packmol_output_file, forcefield, final_box_size,
                                       topologies=topologies, constituent_molecules=constituent_molecules,
                                       mols_number=mols_number, user_settings=setting,
                                       ff_site_property=ff_site_property, input_filename="lammps.in",
                                       data_filename=data_filename, lammps_cmd=lammps_cmd,
                                       db_file=db_file, log_filename=log_filename,
                                       dump_filenames=dump_filenames, parents=[fw_packmol])
        fws.append(fw_lammps)

    return Workflow(fws, name=name)
