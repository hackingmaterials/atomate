# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines workflows for bandstructure calculations.
"""

from fireworks import Workflow
import os
from pymatgen.io.vasp.sets import MPVaspInputSet

from monty.serialization import loadfn
from matmethods.utils.loaders import get_wf_from_spec_dict

__author__ = 'Anubhav Jain, Kiran Mathew'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov'


def get_wf_bandstructure(structure, vasp_input_set=None, vasp_cmd="vasp",
                         db_file=None, add_hse_gap=False):
    """
    Return vasp workflow consisting of 4 fireworks:

    Firework 1 : write vasp input set for structural relaxation,
                 run vasp,
                 pass run location,
                 database insertion.

    Firework 2 : copy files from previous run,
                 write vasp input set for static run,
                 run vasp,
                 pass run location
                 database insertion.

    Firework 3 : copy files from previous run,
                 write vasp input set for non self-consistent
                 (constant charge density) run in
                 uniform mode,
                 run vasp,
                 pass run location
                 database insertion.

    Firework 4 : copy files from previous run,
                 write vasp input set for non self-consistent
                 (constant charge density) run in
                 line mode,
                 run vasp,
                 pass run location
                 database insertion.

    Firework 5 : (optional) HSE gap run

    Args:
        structure (Structure): input structure to be optimized and run
        vasp_input_set (DictVaspInputSet): vasp input set.
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.
        add_hse_gap (bool): add an HSE gap

    Returns:
        Workflow
    """

    wf_file = "band_structure_hsegap.yaml" if add_hse_gap \
        else "band_structure.yaml"

    d = loadfn(os.path.join(os.path.dirname(__file__), wf_file))

    v = vasp_input_set or MPVaspInputSet(force_gamma=True)
    d["fireworks"][0]["params"] = {"vasp_input_set": v.as_dict()}

    d["common_params"] = {
        "vasp_cmd": vasp_cmd,
        "db_file": db_file
    }

    return get_wf_from_spec_dict(structure, d)

if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest

    structure = PymatgenTest.get_structure("Si")
    wf = get_wf_bandstructure(structure)
