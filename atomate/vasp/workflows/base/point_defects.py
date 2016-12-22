# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the point defects workflow: structure optimization followed by transmuter fireworks.
"""

from fireworks import Workflow

from atomate.utils.utils import get_logger
from atomate.vasp.fireworks.core import OptimizeFW, TransmuterFW

from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


def get_wf_point_defects(structure, defect_transformations, defect_transformations_params, charge_states=None,
                         name="point_defects", vasp_input_set=None, lepsilon=False, vasp_cmd="vasp",
                         db_file=None, user_kpoints_settings=None, tag=""):
    """
    Returns a point defects workflow.

    Firework 1 : structural relaxation

    Firework 2 - n_defect_transformations * n_charge_states: Apply the defect transforamtions and
                                    run static calculation for each charge state.


    Args:
        structure (Structure): input structure to be optimized and run
        defect_transformations ([str]): list of defect tranforamtionc lass names as defined in
            pymatgen.transforamtions.defect_transformations.py.
            Example: ["VacancyTransformation", "InterstitialTransformation"]
        defect_transformations_params (list(dict)): list of paramter dictionaries for each defect
            deformation specified in defect_transformations.
            Example: for defect_transformations =  ["VacancyTransformation"],
             defect_transformations_params = [{"supercell_dim":[2,2,2], "species"=None,
                                               "valences"=None, "radii"=None}]
        name (str): some appropriate name for the transmuter fireworks.
        vasp_input_set (DictVaspInputSet): vasp input set for relaxation.
        lepsilon (bool): whether or not compute static dielectric constant/normal modes
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.
        user_kpoints_settings (dict): example: {"grid_density": 7000}
        tag (str): some unique string that will be appended to the names of the fireworks so that
            the data from those tagged fireworks can be queried later during the analysis.

    Returns:
        Workflow
    """
    charge_states = charge_states or []

    # input set for relaxation
    vis_relax = vasp_input_set or MPRelaxSet(structure, force_gamma=True)
    if user_kpoints_settings:
        v = vis_relax.as_dict()
        v.update({"user_kpoints_settings": user_kpoints_settings})
        vis_relax = vis_relax.__class__.from_dict(v)

    # Structure optimization firework
    fws = [OptimizeFW(structure=structure, vasp_input_set=vis_relax, vasp_cmd=vasp_cmd,
                      db_file=db_file, name="{} structure optimization".format(tag))]

    # user input for the static inputset
    uis_static = {"ISTART": 1}
    nelect_input = MPStaticSet(structure).nelect

    # Apply each defect transformation and run static calculations for each charge state
    # ==> n_fireworks = n_defect_transformations * n_charge_states
    for i, dt in enumerate(defect_transformations):
        scell = defect_transformations_params[i]["supercell_dim"]
        ncells = scell[0] * scell[1] * scell[2]
        nelect = nelect_input * defect_transformations_params[i] * ncells
        for cs in charge_states:
            uis_static["NELECT"] = nelect+cs
            vis_static = MPStaticSet(structure, force_gamma=True, lepsilon=lepsilon,
                                     user_kpoints_settings=user_kpoints_settings,
                                     user_incar_settings=uis_static)
            fw = TransmuterFW(name="{} {} {}".format(tag, name, dt), structure=structure,
                              transformations=[dt],
                              transformation_params=[defect_transformations_params[i]],
                              vasp_input_set=vis_static, copy_vasp_outputs=True, parents=fws[0],
                              vasp_cmd=vasp_cmd, db_file=db_file)

            fws.append(fw)

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)

    return Workflow(fws, name=wfname)
