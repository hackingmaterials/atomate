# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the elastic workflow
"""

import json
from decimal import Decimal

import numpy as np

from fireworks import FireTaskBase, Firework, FWAction, Workflow
from fireworks.utilities.fw_serializers import DATETIME_HANDLER
from fireworks.utilities.fw_utilities import explicit_serialize

from matmethods.utils.utils import env_chk, get_logger
from matmethods.vasp.drones import VaspDrone
from matmethods.vasp.database import MMDb
from matmethods.vasp.fireworks.core import OptimizeFW, TransmuterFW

from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.analysis.elasticity.strain import Deformation, IndependentStrain
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.analysis.elasticity import reverse_voigt_map
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.sets import MPRelaxSet, DictSet
from pymatgen import Structure

from matmethods.vasp.fireworks.core import LcalcpolFW
from matmethods.vasp.fireworks.core import HSEBSFW

__author__ = 'Tess Smidt'
__email__ = 'tsmidt@berkeley.edu'

logger = get_logger(__name__)

def get_wf_ferroelectric(polar_structure,nonpolar_structure, pair_id = None, vasp_input_set_polar=None,
                         vasp_input_set_nonpolar=None, vasp_cmd="vasp", db_file=None, nimages = 5):

    """

    Args:
        polar_structure (Structure): polar structure of candidate ferroelectric
        nonpolar_structure (Structure): nonpolar structure of candidate ferroelectric
        vasp_input_set (DictVaspInputSet): vasp input set.
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.
        nimages: Number of interpolations calculated between polar and nonpolar structures.

    Returns:

    """

    # First run polarization calculation on polar structure. Defuse children fireworks if metallic.
    polar = LcalcpolFW(polar_structure,name="polar_polarization",static_name="polar_static", parents=None,
                       vasp_cmd=vasp_cmd,db_file=db_file,vasp_input_set=vasp_input_set_polar,defuse_children=True)


    # Then run polarization calculation on nonpolarstructure structure.
    nonpolar = LcalcpolFW(nonpolar_structure,name="nonpolar_polarization",static_name="nonpolar_static", parents=polar,
                          vasp_cmd=vasp_cmd,db_file=db_file,vasp_input_set=vasp_input_set_nonpolar)


    # Interpolation polarization
    interpolation = []
    for i in range(nimages)[1:-1]:
        interpolation.append(
            LcalcpolFW(nonpolar_structure, name="interpolation_{i}_polarization".format(i=i),
                       static_name="nonpolar_static", vasp_cmd=vasp_cmd, db_file=db_file, vasp_input_set=vasp_input_set,
                       interpolate=True, start="polar_static",end="nonpolar_static",nimages=5,this_image=i,
                       parents=[polar,nonpolar]))

    # Run HSE calculation at band gap for polar calculation if polar structure is not metallic
    hse = HSEBSFW(polar_structure,polar,name="polar_hse_gap",vasp_cmd=vasp_cmd,
                  db_file=db_file,calc_loc="polar_polarization")

    return Workflow([polar,nonpolar]+interpolation+[hse])
