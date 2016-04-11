from matmethods.vasp.input_sets import StructureOptimizationVaspInputSet
from matmethods.vasp.vasp_powerups import use_custodian, decorate_write_name
from matmethods.vasp.workflows.base.band_structure import get_wf_bandstructure_Vasp

__author__ = 'Anubhav Jain <ajain@lbl.gov>'



def wf_band_structure(structure):

    wf = get_wf_bandstructure_Vasp(structure, vasp_input_set=StructureOptimizationVaspInputSet(),
                                   vasp_cmd=">>vasp_cmd<<", db_file=">>db_file<<")

    wf = use_custodian(wf)
    wf = use_custodian(wf, fw_name_constraint="structure optimization",
                       custodian_params={"job_type": "double_relaxation_run"})
    wf = decorate_write_name(wf)

    return wf
