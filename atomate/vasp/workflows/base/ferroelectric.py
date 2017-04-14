# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the elastic workflow
"""

import json
from decimal import Decimal

import numpy as np

from fireworks import FiretaskBase, Firework, FWAction, Workflow

from atomate.utils.utils import env_chk, get_logger

from pymatgen import Structure

from atomate.vasp.fireworks.core import LcalcpolFW, OptimizeFW
from atomate.vasp.fireworks.core import HSEBSFW

__author__ = 'Tess Smidt'
__email__ = 'tsmidt@berkeley.edu'

logger = get_logger(__name__)


def get_wf_ferroelectric(polar_structure,nonpolar_structure, vasp_cmd="vasp", db_file=None,
                         vasp_input_set_polar=None, vasp_input_set_nonpolar=None,
                         relax=False, vasp_relax_input_set_polar=None, vasp_relax_input_set_nonpolar=None,
                         nimages = 5, hse = False, task_prefix = "",from_prev_settings=None):
    """
    Returns a workflow to calculate the spontaneous polarization of polar_structure using
    a nonpolar reference phase structure and linear interpolations between the polar and
    nonpolar structure.

    Args:
        polar_structure (Structure): polar structure of candidate ferroelectric
        nonpolar_structure (Structure): nonpolar structure of candidate ferroelectric
        vasp_input_set_polar (DictVaspInputSet): vasp polar input set.
        vasp_input_set_nonpolar (DictVaspInputSet): vasp nonpolar input set.
        vasp_relax_input_set_polar (DictVaspInputSet): vasp polar input set.
        vasp_relax_input_set_nonpolar (DictVaspInputSet): vasp nonpolar input set.
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.
        nimages: Number of interpolations calculated between polar and nonpolar structures.

    Returns:

    """
    wf = []

    if relax:
        polar_relax = OptimizeFW(polar_structure,name="{t}polar_relaxation".format(t=task_prefix),
                                 vasp_cmd=vasp_cmd, db_file=db_file, vasp_input_set=vasp_relax_input_set_polar)
        nonpolar_relax = OptimizeFW(nonpolar_structure, name="{t}nonpolar_relaxation".format(t=task_prefix),
                                    vasp_cmd=vasp_cmd, db_file=db_file, vasp_input_set=vasp_relax_input_set_nonpolar)
        wf.append(polar_relax)
        wf.append(nonpolar_relax)
        parents_polar = polar_relax
        parents_nonpolar = nonpolar_relax
        write_from_prev = True
        calc_loc_polar = "{t}polar_relaxation".format(t=task_prefix)
        calc_loc_nonpolar = "{t}nonpolar_relaxation".format(t=task_prefix)
    else:
        parents_polar = None
        parents_nonpolar = None
        write_from_prev = False
        calc_loc_polar = True
        calc_loc_nonpolar = True

    # First run polarization calculation on polar structure. Defuse children fireworks if metallic.
    polar = LcalcpolFW(polar_structure,
                       name="{t}polar_polarization".format(t=task_prefix),
                       static_name="{t}polar_static".format(t=task_prefix),
                       parents=parents_polar,
                       vasp_cmd=vasp_cmd,db_file=db_file,
                       vasp_input_set=vasp_input_set_polar,
                       write_from_prev=write_from_prev,
                       defuse_children=True,calc_loc=calc_loc_polar,from_prev_settings=from_prev_settings)


    # Then run polarization calculation on nonpolarstructure structure.
    nonpolar = LcalcpolFW(nonpolar_structure,
                          name="{t}nonpolar_polarization".format(t=task_prefix),
                          static_name="{t}nonpolar_static".format(t=task_prefix),
                          parents=parents_nonpolar,
                          vasp_cmd=vasp_cmd,db_file=db_file,
                          vasp_input_set=vasp_input_set_nonpolar,
                          write_from_prev=write_from_prev,defuse_children=True,
                          calc_loc=calc_loc_nonpolar,from_prev_settings=from_prev_settings)

    # Interpolation polarization
    interpolation = []
    # this is to start from one increment after polar and end prior to nonpolar
    # the Structure.interpolate method adds an additional image for the nonpolar endpoint
    for i in range(1,nimages):
        interpolation.append(
            LcalcpolFW(nonpolar_structure,
                       name="{t}interpolation_{i}_polarization".format(i=i,t=task_prefix),
                       static_name="{t}interpolation_{i}_static".format(i=i,t=task_prefix),
                       vasp_cmd=vasp_cmd, db_file=db_file,
                       vasp_input_set=vasp_input_set_polar,
                       write_from_prev=False, interpolate=True,
                       start="{t}polar_static".format(t=task_prefix),
                       end="{t}nonpolar_static".format(t=task_prefix),
                       nimages=nimages,this_image=i, parents=[polar,nonpolar],
                       from_prev_settings=from_prev_settings))

    wf.append(polar)
    wf.append(nonpolar)
    wf += interpolation

    # Add FireTask that uses Polarization object to store spontaneous polarization information

    # Delete?
    if hse:
        # Run HSE calculation at band gap for polar calculation if polar structure is not metallic
        hse = HSEBSFW(polar_structure,polar,name="{t}polar_hse_gap".format(t=task_prefix),vasp_cmd=vasp_cmd,
                      db_file=db_file,calc_loc="{t}polar_polarization".format(t=task_prefix))
        wf.append(hse)

    return Workflow(wf)
