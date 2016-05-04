from pymatgen import Structure
from pymatgen.io.vasp.sets import MPVaspInputSet

from fireworks import LaunchPad

from matmethods.vasp.workflows.base.spinorbit_coupling import get_wf_spinorbit_coupling


if __name__ == "__main__":
    # the structure from vasp wiki example
    fe_monomer = Structure([[1.73, 1.73, 0.0],
                            [-1.73, 1.73, 0.0],
                            [0.0, 0.0, 10.0]],
                           ["Fe"],
                           [[0, 0, 0]])
    user_incar_settings = {"ISIF": 2,
                           "NPAR": 4,
                           "ALGO": "Normal",
                           "LREAL": ".FALSE."}
    vis = MPVaspInputSet(user_incar_settings=user_incar_settings, force_gamma=True)
    wf = get_wf_spinorbit_coupling(fe_monomer, [3.0], field_directions=[[0,0,1]],
                                   vasp_input_set=vis, vasp_cmd="srun vasp",
                                   vasp_ncl="srun vasp_ncl",  db_file=">>db_file<<")
