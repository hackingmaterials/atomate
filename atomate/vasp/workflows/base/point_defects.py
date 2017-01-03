# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the point defects workflow: structure optimization followed by transmuter fireworks.
"""

from fireworks import Workflow

from atomate.utils.utils import get_logger
from atomate.vasp.fireworks.core import OptimizeFW, TransmuterFW, StaticFW

from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


def get_wf_point_defects(structure, defect_transformations, defect_transformations_params,
                         charge_states=None, name="point_defects", vasp_input_set=None,
                         lepsilon=False, vasp_cmd="vasp", db_file=None, user_kpoints_settings=None,
                         tag=""):
    """
    Returns a point defects workflow.

    Firework 1 : structural relaxation

    Firework 2 : static

    Firework 3 - : Apply the defect transformations and run static calculation for each charge state.


    Args:
        structure (Structure): input structure to be optimized and run
        defect_transformations ([str]): list of defect transformation class names as defined in
            pymatgen.transforamtions.defect_transformations.py.
            Example: ["VacancyTransformation", "InterstitialTransformation"]
        defect_transformations_params (list(dict)): list of parameter dictionaries for each defect
            deformation specified in defect_transformations.
            Example: for defect_transformations =  ["VacancyTransformation"],
             defect_transformations_params = [{"supercell_dim":[2,2,2], "species":None,
                                               "valences":None, "radii":None}]
        charge_states (list): list of floats indicating the extra number of electrons to be added
            to the system(updates NELECT).
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

    TODO: insert task to compute average electrostatic potential and pass it along
    """
    # if not given, just do the neutral system calculations
    charge_states = charge_states or [0]

    # input set for relaxation
    vis_relax = vasp_input_set or MPRelaxSet(structure, force_gamma=True)
    if user_kpoints_settings:
        v = vis_relax.as_dict()
        v.update({"user_kpoints_settings": user_kpoints_settings})
        vis_relax = vis_relax.__class__.from_dict(v)

    # Structure optimization firework
    fws = [OptimizeFW(structure=structure, vasp_input_set=vis_relax, vasp_cmd=vasp_cmd,
                      db_file=db_file, name="{} structure optimization".format(tag))]

    # Static firework: copy previous outputs and do a static run on the defect-free system
    uis_static = {"ISTART": 1}
    uis_static_relax_ions = {"ISTART": 1, "IBRION": 2, "NSW": 100}
    vis_static = MPStaticSet(structure, force_gamma=True, lepsilon=lepsilon,
                             user_kpoints_settings=user_kpoints_settings,
                             user_incar_settings=uis_static_relax_ions)
    nelect_input = vis_static.nelect
    fws.append(StaticFW(structure, vasp_input_set=vis_static, vasp_cmd=vasp_cmd,
                        db_file=db_file, name="{} static calculation".format(tag),
                        copy_vasp_outputs=True, parents=fws[0]))

    # Apply each defect transformation and run static calculations for each charge state
    # ==> n_fireworks = n_defect_transformations * n_charge_states + n_defect_transformations
    for i, dt in enumerate(defect_transformations):
        # Transmuter firework: copy previous outputs, create supercell with the defect(neutral) and
        # relax the ions
        vis_static = MPStaticSet(structure, force_gamma=True, lepsilon=lepsilon,
                                 user_kpoints_settings=user_kpoints_settings,
                                 user_incar_settings=uis_static_relax_ions)
        fw_neutral_defect = TransmuterFW(name="{}_{}_{}".format(tag, name, dt),  structure=structure,
                                transformations=[dt],
                                transformation_params=[defect_transformations_params[i]],
                                vasp_input_set=vis_static, copy_vasp_outputs=True, parents=fws[1],
                                vasp_cmd=vasp_cmd, db_file=db_file)
        fws.append(fw_neutral_defect)

        # compute the nelect for the supercell
        scell = defect_transformations_params[i]["supercell_dim"]
        ncells = scell[0] * scell[1] * scell[2]
        nelect = nelect_input * ncells

        for cs in charge_states:
            uis_static["NELECT"] = nelect+cs
            vis_static = MPStaticSet(structure, force_gamma=True, lepsilon=lepsilon,
                                     user_kpoints_settings=user_kpoints_settings,
                                     user_incar_settings=uis_static)
            # Static firework: copy outputs from neutral defect calculation and do a static
            # calculation for the given charge state
            fws.append(StaticFW(structure, vasp_input_set=vis_static, vasp_cmd=vasp_cmd,
                                db_file=db_file, parents=fw_neutral_defect,
                                name="{} {} {} static calculation".format(tag, dt, cs)))

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)

    return Workflow(fws, name=wfname)
