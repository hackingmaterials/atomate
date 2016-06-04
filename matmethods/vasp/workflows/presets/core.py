from matmethods.vasp.workflows.base.optical import get_wf_dielectric_constant

from matmethods.vasp.workflows.base.core import get_wf
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet

from matmethods.vasp.vasp_powerups import add_namefile, \
    add_small_gap_multiply, use_scratch_dir

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

# TODO: allow standard behavior to be modified via a config file (but retain structure as only param)
# TODO: add SOC workflow


def wf_band_structure(structure):
    wf = get_wf(structure, "band_structure.yaml", vasp_cmd=">>vasp_cmd<<",
                db_file=">>db_file<<")
    wf = add_namefile(wf)
    wf = add_small_gap_multiply(wf, 0.5, 5, "static")
    wf = add_small_gap_multiply(wf, 0.5, 5, "nscf")
    wf = use_scratch_dir(wf, ">>scratch_dir<<")
    return wf


def wf_band_structure_plus_hse(structure):
    wf = get_wf(structure, "band_structure_hsegap.yaml", vasp_cmd=">>vasp_cmd<<",
                db_file=">>db_file<<")
    wf = add_namefile(wf)
    wf = add_small_gap_multiply(wf, 0.5, 5, "static")
    wf = add_small_gap_multiply(wf, 0.5, 5, "nscf")
    wf = use_scratch_dir(wf, ">>scratch_dir<<")
    return wf


def wf_static(structure):
    wf = get_wf(structure, "static_only.yaml", vasp_cmd=">>vasp_cmd<<",
                db_file=">>db_file<<")

    wf = add_namefile(wf)
    wf = use_scratch_dir(wf, ">>scratch_dir<<")

    return wf


def wf_structure_optimization(structure):
    wf = get_wf(structure, "optimize_only.yaml", vasp_cmd=">>vasp_cmd<<",
                db_file=">>db_file<<")

    wf = add_namefile(wf)
    wf = use_scratch_dir(wf, ">>scratch_dir<<")

    return wf


def wf_dielectric_constant(structure):
    wf = get_wf(structure, "dielectric_constant.yaml", vasp_cmd=">>vasp_cmd<<",
                db_file=">>db_file<<")

    wf = add_namefile(wf)
    wf = use_scratch_dir(wf, ">>scratch_dir<<")

    return wf