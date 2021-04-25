__author__ = "Anubhav Jain <ajain@lbl.gov>"

# TODO: @computron should be able to load from YAML -computron

ADD_NAMEFILE = True
SCRATCH_DIR = ">>scratch_dir<<"
GAMMA_VASP_CMD = ">>gamma_vasp_cmd<<"
SMALLGAP_KPOINT_MULTIPLY = True
ADD_MODIFY_INCAR = False
STABILITY_CHECK = False
VASP_CMD = ">>vasp_cmd<<"
VDW_KERNEL_DIR = ">>vdw_kernel_dir<<"
DB_FILE = ">>db_file<<"
ADD_WF_METADATA = True
LOBSTER_CMD = ">>lobster_cmd<<"

# whether to use only half the kpoint density in
# the initial relaxation of a structure optimization for faster performance
HALF_KPOINTS_FIRST_RELAX = False

# maximum force allowed on atom for successful structure optimization
RELAX_MAX_FORCE = 0.25

# this is a three-way toggle on what to do if your job looks OK,
# but is actually unconverged (either electronic or ionic).
# True -> mark job as COMPLETED, but defuse children.
# False --> do nothing, continue with workflow as normal.
# "fizzle" --> throw an error (mark this job as FIZZLED)
DEFUSE_UNSUCCESSFUL = "fizzle"

# maximum number of errors to correct before custodian gives up
CUSTODIAN_MAX_ERRORS = 5

# store data from these files in database if present
STORE_VOLUMETRIC_DATA = ()  # e.g. ("chgcar", "aeccar0", "aeccar2", "elfcar", "locpot")

# ingest any additional JSON data present into database when parsing VASP directories
# useful for storing duplicate of FW.json
STORE_ADDITIONAL_JSON = False

# vasp output files that will be copied to lobster run
VASP_OUTPUT_FILES = [
    "OUTCAR",
    "vasprun.xml",
    "CHG",
    "CHGCAR",
    "CONTCAR",
    "INCAR",
    "KPOINTS",
    "POSCAR",
    "POTCAR",
    "DOSCAR",
    "EIGENVAL",
    "IBZKPT",
    "OSZICAR",
    "WAVECAR",
    "XDATCAR",
]

# lobster output files
LOBSTERINPUT_FILES = ["lobsterin"]
LOBSTEROUTPUT_FILES = [
    "lobsterout",
    "CHARGE.lobster",
    "COHPCAR.lobster",
    "COOPCAR.lobster",
    "DOSCAR.lobster",
    "GROSSPOP.lobster",
    "ICOHPLIST.lobster",
    "ICOOPLIST.lobster",
    "lobster.out",
    "projectionData.lobster",
]
