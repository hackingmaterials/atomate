__author__ = 'Anubhav Jain <ajain@lbl.gov>'

# TODO: @computron should be able to load from YAML -computron

ADD_NAMEFILE = True
SCRATCH_DIR = ">>scratch_dir<<"
GAMMA_VASP_CMD = ">>gamma_vasp_cmd<<"
SMALLGAP_KPOINT_MULTIPLY = True
ADD_MODIFY_INCAR = False
STABILITY_CHECK = False
VASP_CMD = ">>vasp_cmd<<"
DB_FILE = ">>db_file<<"
ADD_WF_METADATA = True
HALF_KPOINTS_FIRST_RELAX = False  # whether to use only half the kpoint density in the initial relaxation of a structure optimization for faster performance
RELAX_MAX_FORCE = 0.25  # maximum force allowed on atom for successful structure optimization