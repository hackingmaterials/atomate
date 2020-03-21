from pymatgen.io.cp2k.sets import (
    StaticSet,
    RelaxSet,
    HybridStaticSet,
    HybridRelaxSet,
)
from atomate.cp2k.fireworks.core import (
    StaticFW,
    RelaxFW,
    StaticHybridFW,
    RelaxHybridFW,
)
from fireworks import Workflow

ADD_NAMEFILE = True
SCRATCH_DIR = ">>scratch_dir<<"
STABILITY_CHECK = False
CP2K_CMD = ">>cp2k_cmd<<"
DB_FILE = ">>db_file<<"
ADD_WF_METADATA = True


def get_wf_static(
    structure,
    cp2k_input_set=None,
    name="Static-WF",
    cp2k_cmd=">>cp2k_cmd<<",
    db_file=">>db_file<<",
    user_cp2k_settings=None,
    metadata=None,
):
    """
    Returns the workflow that computes the bulk modulus by fitting to the given equation of state.

    Args:
        structure (Structure): input structure.
        cp2k_input_set (Cp2kInputSet): Cp2k input set to use, if not using the default.
        cp2k_cmd (str): cp2k command to run.
        db_file (str): path to the db file.
        tag (str): something unique to identify the tasks in this workflow. If None a random uuid
            will be assigned.
        user_cp2k_settings (dict): If passing a non-default Cp2kInputSet, this dict contains the
            kwargs to pass to the Cp2kInputSet initialization. Note that to override the cp2k input
            file parameters themselves, a key of the form {'override_default_params': {...}}

    Returns:
        Workflow
    """
    fws = []

    cis_static = cp2k_input_set or StaticSet(structure, **user_cp2k_settings)

    fw = StaticFW(
        structure=structure,
        name=name,
        cp2k_input_set=cis_static,
        cp2k_input_set_params=user_cp2k_settings,
        cp2k_cmd=cp2k_cmd,
        prev_calc_loc=False,
        prev_calc_dir=None,
        db_file=db_file,
        cp2ktodb_kwargs=None,
        parents=None,
    )
    fws.append(fw)

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)

    return Workflow(fws, name=wfname, metadata=metadata)


def get_wf_relax(
    structure,
    cp2k_input_set=None,
    name="Relax WF",
    cp2k_cmd=CP2K_CMD,
    db_file=DB_FILE,
    user_cp2k_settings=None,
    metadata=None,
):
    """
    Returns the workflow that computes the bulk modulus by fitting to the given equation of state.

    Args:
        structure (Structure): input structure.
        cp2k_input_set (Cp2kInputSet): for relax calculations
        cp2k_cmd (str): cp2k command to run.
        db_file (str): path to the db file.
        tag (str): something unique to identify the tasks in this workflow. If None a random uuid
            will be assigned.

    Returns:
        Workflow
    """
    fws = []

    cis_static = cp2k_input_set or RelaxSet(structure, **user_cp2k_settings)

    fw = RelaxFW(
        structure=structure,
        name=name,
        cp2k_input_set=cis_static,
        cp2k_input_set_params=user_cp2k_settings,
        cp2k_cmd=cp2k_cmd,
        prev_calc_loc=False,
        prev_calc_dir=None,
        db_file=db_file,
        cp2ktodb_kwargs=None,
        parents=None,
    )
    fws.append(fw)

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)

    return Workflow(fws, name=wfname, metadata=metadata)


def get_wf_hybrid_static(
    structure,
    cp2k_static_input_set=None,
    cp2k_hybrid_input_set=None,
    name="Hybrid-Static-WF",
    user_static_settings={},
    user_hybrid_settings={},
    cp2k_cmd=CP2K_CMD,
    db_file=DB_FILE,
    cp2ktodb_kwargs=None,
    metadata=None,
):

    fws = []

    # TODO I really don't like this work around... currently I'm asserting that all cp2k input files
    # must have the same project name, that way its easier for different fws to find the files from
    # previous fireworks. Should be more flexible. -NW
    if "project_name" not in user_static_settings.keys():
        user_static_settings["project_name"] = name
    if "project_name" not in user_hybrid_settings.keys():
        user_hybrid_settings["project_name"] = name

    cp2k_static_input_set = cp2k_static_input_set or StaticSet(
        structure, **user_static_settings
    )
    cp2k_hybrid_input_set = cp2k_hybrid_input_set or HybridStaticSet(
        structure, **user_hybrid_settings
    )

    fw1 = StaticFW(
        structure=structure,
        name=name,
        cp2k_input_set=cp2k_static_input_set,
        cp2k_cmd=cp2k_cmd,
        prev_calc_loc=False,
        db_file=db_file,
        cp2ktodb_kwargs=cp2ktodb_kwargs,
        parents=None,
    )
    fws.append(fw1)

    fw2 = StaticHybridFW(
        structure=structure,
        name=name,
        cp2k_input_set=cp2k_hybrid_input_set,
        cp2k_cmd=cp2k_cmd,
        prev_calc_loc=True,
        db_file=db_file,
        cp2ktodb_kwargs=cp2ktodb_kwargs,
        parents=fw1,
    )
    fws.append(fw2)

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)

    return Workflow(fws, name=wfname, metadata=metadata)
