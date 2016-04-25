from pymatgen import Structure

from fireworks import LaunchPad

from matmethods.vasp.workflows.base.spinorbit_coupling import get_wf_spinorbit_coupling
from matmethods.vasp.input_sets import StructureOptimizationVaspInputSet


if __name__ == "__main__":
    
    fe_monomer = Structure([[1.73, 1.73,0.0],[-1.73,1.73,0.0],[0.0,0.0,8.63]], 
                         ["Fe"], [[0,0,0]])
    config_dict_override = {"INCAR":{"ISIF":2}}
    vis = StructureOptimizationVaspInputSet(config_dict_override=config_dict_override)
    lp = LaunchPad.auto_load()
    wf = get_wf_spinorbit_coupling(fe_monomer, [3.0], field_directions=[[0,0,1]],
                                   vasp_input_set = vis,
                                   vasp_cmd="srun -n 16 vasp", vasp_ncl="srun -n 16 vasp_ncl",  
                                   db_file=">>db_file<<")
    lp.add_wf(wf)
