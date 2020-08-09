from pymatgen.io.cp2k.sets import (
    StaticSet,
    RelaxSet,
    HybridStaticSet,
    HybridRelaxSet,
    HybridCellOptSet
)
from atomate.cp2k.fireworks.core import (
    StaticFW,
    RelaxFW,
    StaticHybridFW,
    RelaxHybridFW,
    CellOptHybridFW
)
from fireworks import Workflow
from atomate.utils.utils import get_wf_from_spec_dict
from monty.serialization import loadfn
import os

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

ADD_NAMEFILE = True
SCRATCH_DIR = ">>scratch_dir<<"
STABILITY_CHECK = False
CP2K_CMD = ">>cp2k_cmd<<"
DB_FILE = ">>db_file<<"
ADD_WF_METADATA = True

# TODO Incomplete. Just an idea.
def get_wf(structure, wf_filename, params=None, common_params=None, wf_metadata=None):
    """
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
    """

    fws = []

    sets = [

    ]

    for p in params:
        fw = p['fw']

        mod = __import__(modulepath, globals(), locals(), [classname], 0)
        return getattr(mod, classname)

    cp2k_static_input_set = cp2k_static_input_set or StaticSet(
        structure, **user_static_settings
    )
    cp2k_hybrid_input_set = cp2k_hybrid_input_set or HybridStaticSet(
        structure, **user_hybrid_settings
    )

    gga_name = "{}-GGA-FW".format(structure.composition.reduced_formula)
    hybrid_name = "{}-Hybrid-FW".format(structure.composition.reduced_formula)
    restart_filename = "{}-RESTART.wfn".format(gga_name)  # GGA restart WFN

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)

    return Workflow(fws, name=wfname, metadata=wf_metadata)


def get_wf_static(
    structure,
    cp2k_input_set=None,
    name="Static-WF",
    cp2k_cmd=">>cp2k_cmd<<",
    db_file=">>db_file<<",
    user_cp2k_settings={},
    cp2ktodb_kwargs=None,
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

    fw = StaticFW(
        structure=structure,
        name=name,
        cp2k_input_set=cp2k_input_set,
        cp2k_input_set_params=user_cp2k_settings,
        cp2k_cmd=cp2k_cmd,
        prev_calc_loc=False,
        prev_calc_dir=None,
        db_file=db_file,
        cp2ktodb_kwargs=cp2ktodb_kwargs,
        parents=None,
        files_to_copy=None
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
    user_cp2k_settings={},
    cp2ktodb_kwargs=None,
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

    fw = RelaxFW(
        structure=structure,
        name=name,
        cp2k_input_set=cp2k_input_set,
        cp2k_input_set_params=user_cp2k_settings,
        cp2k_cmd=cp2k_cmd,
        prev_calc_loc=False,
        prev_calc_dir=None,
        db_file=db_file,
        cp2ktodb_kwargs=cp2ktodb_kwargs,
        parents=None,
        files_to_copy=None
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

    gga_name = "{}-GGA-FW".format(structure.composition.reduced_formula)
    hybrid_name = "{}-Hybrid-FW".format(structure.composition.reduced_formula)
    restart_filename = "{}-RESTART.wfn".format(gga_name)  # GGA restart WFN
    user_hybrid_settings['wfn_restart_file_name'] = restart_filename

    fw1 = StaticFW(
        structure=structure,
        name=gga_name,
        cp2k_input_set=cp2k_static_input_set,
        cp2k_input_set_params=user_static_settings,
        cp2k_cmd=cp2k_cmd,
        prev_calc_loc=False,
        db_file=db_file,
        cp2ktodb_kwargs=cp2ktodb_kwargs,
        parents=None,
        files_to_copy=None
    )
    fws.append(fw1)

    fw2 = StaticHybridFW(
        structure=structure,
        name=hybrid_name,
        cp2k_input_set=cp2k_hybrid_input_set,
        cp2k_input_set_params=user_hybrid_settings,
        cp2k_cmd=cp2k_cmd,
        prev_calc_loc=True,
        db_file=db_file,
        cp2ktodb_kwargs=cp2ktodb_kwargs,
        parents=fw1,
        files_to_copy=restart_filename
    )
    fws.append(fw2)

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)

    return Workflow(fws, name=wfname, metadata=metadata)


def get_wf_hybrid_relax(
    structure,
    cp2k_static_input_set=None,
    cp2k_hybrid_input_set=None,
    name="Hybrid-Relax-WF",
    user_static_settings={},
    user_hybrid_settings={},
    cp2k_cmd=CP2K_CMD,
    db_file=DB_FILE,
    cp2ktodb_kwargs=None,
    metadata=None,
):

    fws = []

    gga_name = "{}-GGA-FW".format(structure.composition.reduced_formula)
    hybrid_name = "{}-Hybrid-FW".format(structure.composition.reduced_formula)
    restart_filename = "{}-RESTART.wfn".format(gga_name)  # GGA restart WFN
    user_hybrid_settings['wfn_restart_file_name'] = restart_filename

    fw1 = StaticFW(
        structure=structure,
        name=gga_name,
        cp2k_input_set=cp2k_static_input_set,
        cp2k_input_set_params=user_static_settings,
        cp2k_cmd=cp2k_cmd,
        prev_calc_loc=False,
        db_file=db_file,
        cp2ktodb_kwargs=cp2ktodb_kwargs,
        parents=None,
        files_to_copy=False
    )
    fws.append(fw1)

    fw2 = RelaxHybridFW(
        structure=structure,
        name=hybrid_name,
        cp2k_input_set=cp2k_hybrid_input_set,
        cp2k_input_set_params=user_hybrid_settings,
        cp2k_cmd=cp2k_cmd,
        prev_calc_loc=True,
        db_file=db_file,
        cp2ktodb_kwargs=cp2ktodb_kwargs,
        parents=fw1,
        files_to_copy=restart_filename
    )
    fws.append(fw2)

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)

    return Workflow(fws, name=wfname, metadata=metadata)


def get_wf_hybrid_cell_opt(
    structure,
    cp2k_static_input_set=None,
    cp2k_hybrid_input_set=None,
    name="Hybrid-CellOpt-WF",
    user_static_settings={},
    user_hybrid_settings={},
    cp2k_cmd=CP2K_CMD,
    db_file=DB_FILE,
    cp2ktodb_kwargs=None,
    metadata=None,
):

    fws = []

    gga_name = "{}-GGA-FW".format(structure.composition.reduced_formula)
    hybrid_name = "{}-Hybrid-FW".format(structure.composition.reduced_formula)
    restart_filename = "{}-RESTART.wfn".format(gga_name)  # GGA restart WFN
    user_hybrid_settings['wfn_restart_file_name'] = restart_filename

    fw1 = StaticFW(
        structure=structure,
        name=gga_name,
        cp2k_input_set=cp2k_static_input_set,
        cp2k_input_set_params=user_static_settings,
        cp2k_cmd=cp2k_cmd,
        prev_calc_loc=False,
        db_file=db_file,
        cp2ktodb_kwargs=cp2ktodb_kwargs,
        parents=None,
        files_to_copy=False
    )
    fws.append(fw1)

    fw2 = CellOptHybridFW(
        structure=structure,
        name=hybrid_name,
        cp2k_input_set=cp2k_hybrid_input_set,
        cp2k_input_set_params=user_hybrid_settings,
        cp2k_cmd=cp2k_cmd,
        prev_calc_loc=True,
        db_file=db_file,
        cp2ktodb_kwargs=cp2ktodb_kwargs,
        parents=fw1,
        files_to_copy=restart_filename
    )
    fws.append(fw2)

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)

    return Workflow(fws, name=wfname, metadata=metadata)


def get_wf_gga_relax_to_hybrid_static(
    structure,
    cp2k_gga_input_set=None,
    cp2k_hybrid_input_set=None,
    name="Hybrid-Static-WF",
    user_gga_settings={},
    user_hybrid_settings={},
    cp2k_cmd=CP2K_CMD,
    db_file=DB_FILE,
    cp2ktodb_kwargs=None,
    metadata=None,
):

    fws = []

    gga_name = "{}-GGA-FW".format(structure.composition.reduced_formula)
    hybrid_name = "{}-Hybrid-FW".format(structure.composition.reduced_formula)
    restart_filename = "{}-RESTART.wfn".format(gga_name)  # GGA restart WFN
    user_hybrid_settings['wfn_restart_file_name'] = restart_filename

    fw1 = RelaxFW(
        structure=structure,
        name=gga_name,
        cp2k_input_set=cp2k_gga_input_set,
        cp2k_input_set_params=user_gga_settings,
        cp2k_cmd=cp2k_cmd,
        prev_calc_loc=False,
        db_file=db_file,
        cp2ktodb_kwargs=cp2ktodb_kwargs,
        parents=None,
        files_to_copy=None
    )
    fws.append(fw1)

    fw2 = StaticHybridFW(
        structure=structure,
        name=hybrid_name,
        cp2k_input_set=cp2k_hybrid_input_set,
        cp2k_input_set_params=user_hybrid_settings,
        cp2k_cmd=cp2k_cmd,
        prev_calc_loc=True,
        db_file=db_file,
        cp2ktodb_kwargs=cp2ktodb_kwargs,
        parents=fw1,
        files_to_copy=restart_filename
    )
    fws.append(fw2)

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)

    return Workflow(fws, name=wfname, metadata=metadata)

def get_wf_atom(
    structure,
    cp2k_gga_input_set=None,
    cp2k_hybrid_input_set=None,
    name="Hybrid-Static-WF",
    user_gga_settings={},
    user_hybrid_settings={},
    cp2k_cmd=CP2K_CMD,
    db_file=DB_FILE,
    cp2ktodb_kwargs=None,
    metadata=None,
):

    fws = []

    hybrid_name = "{}-Hybrid-Atom-FW".format(structure.composition.reduced_formula)

    cis = HybridStaticSet(structure, hybrid_functional='PBE0', hf_fraction=1, gga_x_fraction=0)
    cis.insert({'FORCE_EVAL': {'DFT': {''}}})

    fw = StaticHybridFW(
        structure=structure,
        name=hybrid_name,
        cp2k_input_set=cp2k_hybrid_input_set,
        cp2k_input_set_params=user_hybrid_settings,
        cp2k_cmd=cp2k_cmd,
        prev_calc_loc=True,
        db_file=db_file,
        cp2ktodb_kwargs=cp2ktodb_kwargs,
        parents=None,
        files_to_copy=None
    )
    fws.append(fw)

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)

    return Workflow(fws, name=wfname, metadata=metadata)
