from matmethods.vasp.input_sets import StructureOptimizationVaspInputSet, StaticVaspInputSet
from matmethods.vasp.vasp_powerups import use_custodian, decorate_write_name, \
    add_small_gap_multiply, use_scratch_dir
from matmethods.vasp.workflows.base.band_structure import get_wf_bandstructure
from matmethods.vasp.workflows.base.single_vasp import get_wf_single

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

# TODO: allow standard behavior to be modified via a config file (but keep structure as only param)


def wf_band_structure(structure):

    """
    optimizes structure, then computes both uniform and line mode band structures
    :param structure:
    :return:
    """
    wf = get_wf_bandstructure(structure, vasp_input_set=StructureOptimizationVaspInputSet(),
                              vasp_cmd=">>vasp_cmd<<", db_file=">>db_file<<")
    wf = use_custodian(wf)
    wf = use_custodian(wf, fw_name_constraint="structure optimization",
                       custodian_params={"job_type": "double_relaxation_run",
                                         "max_force_threshold": 0.25})
    wf = decorate_write_name(wf)
    wf = add_small_gap_multiply(wf, 0.5, 5, "static")
    wf = add_small_gap_multiply(wf, 0.5, 5, "nscf")
    wf = use_scratch_dir(wf, ">>scratch_dir<<")

    return wf


def wf_static(structure):
    """
    single static calculation
    :param structure:
    :return:
    """
    wf = get_wf_single(structure, vasp_input_set=StaticVaspInputSet(), vasp_cmd=">>vasp_cmd<<",
                       db_file=">>db_file<<", task_label="static")
    wf = use_custodian(wf)
    wf = decorate_write_name(wf)
    wf = use_scratch_dir(wf, ">>scratch_dir<<")

    return wf


def wf_structure_optimization(structure):
    """
    single structure optimization calculation
    :param structure:
    :return:
    """
    wf = get_wf_single(structure, vasp_input_set=StructureOptimizationVaspInputSet(),
                       vasp_cmd=">>vasp_cmd<<", db_file=">>db_file<<",
                       task_label="structure optimization")
    wf = use_custodian(wf, fw_name_constraint="structure optimization",
                       custodian_params={"job_type": "double_relaxation_run",
                                         "max_force_threshold": 0.25})
    wf = decorate_write_name(wf)
    wf = use_scratch_dir(wf, ">>scratch_dir<<")

    return wf
