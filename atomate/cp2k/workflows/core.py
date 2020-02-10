from pymatgen.io.cp2k.sets import StaticSet
from atomate.cp2k.fireworks.core import StaticFW
from fireworks import Workflow

ADD_NAMEFILE = True
SCRATCH_DIR = ">>scratch_dir<<"
GAMMA_VASP_CMD = ">>gamma_vasp_cmd<<"
SMALLGAP_KPOINT_MULTIPLY = True
ADD_MODIFY_INCAR = False
STABILITY_CHECK = False
VASP_CMD = ">>cp2k_cmd<<"
DB_FILE = ">>db_file<<"
ADD_WF_METADATA = True


def get_wf_static(structure, cp2k_input_set=None, name='Static WF',
                  cp2k_cmd="cp2k.popt", db_file=None,
                  user_cp2k_settings=None, metadata=None):
    """
    Returns the workflow that computes the bulk modulus by fitting to the given equation of state.

    Args:
        structure (Structure): input structure.
        deformations (list): list of deformation matrices(list of lists).
        vasp_input_set (VaspInputSet): for the static deformation calculations
        vasp_cmd (str): vasp command to run.
        db_file (str): path to the db file.
        user_kpoints_settings (dict): example: {"grid_density": 7000}
        eos (str): equation of state used for fitting the energies and the volumes.
            supported equation of states: "quadratic", "murnaghan", "birch", "birch_murnaghan",
            "pourier_tarantola", "vinet", "deltafactor". See pymatgen.analysis.eos.py
        tag (str): something unique to identify the tasks in this workflow. If None a random uuid
            will be assigned.
        user_incar_settings (dict):

    Returns:
        Workflow
    """
    fws = []

    cis_static = cp2k_input_set or StaticSet(structure)

    fw = StaticFW(structure=structure, name=name, cp2k_input_set=cis_static,
                  cp2k_input_set_params=user_cp2k_settings, cp2k_cmd=cp2k_cmd,
                  prev_calc_loc=True, prev_calc_dir=None, db_file=None,
                  cp2ktodb_kwargs=None, parents=None)
    fws.append(fw)

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)

    return Workflow(fws, name=wfname, metadata=metadata)

