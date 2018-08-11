# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

from fireworks.core.firework import Workflow
from atomate.vasp.fireworks.core import OptimizeFW, TransmuterFW
from atomate.vasp.powerups import add_additional_fields_to_taskdocs, add_tags
from pymatgen.io.vasp.sets import MVLGBSet

__author__ = 'Hui Zheng'
__email__ = 'huz071@eng.ucsd.edu'


def get_gb_fw(bulk_structure, gb_gen_params=None, db_file=None,
              vasp_input_set=None, parents=None, vasp_cmd="vasp", name=None):
    """
    The firework will generate a grain boundary from given bulk structure and gb_gen_params
     and then return a TransmuterFW which include the relaxation for generated GB structure.

    Args:
        bulk_structure (Structure): Bulk structure could be primitive cell or conventional unit cell.
        gb_gen_params (dict): dictionary of gb generation parameters, used to generate gb from
            GBGenerator (which could be called by: from pymatgen.analysis.gb.gb import GBGenerator,
            where there are details of description for each parameter)
            gb_gen_params is necessary to get the gb that corresponds to the bulk structure.
            One example is given below:
              gb_s5_001_params = {"rotation_axis": [1, 0, 0], "rotation_angle": 36.86989764584402,
                        "expand_times": 2, "vacuum_thickness": 0.0, "normal": True,
                        "ratio": None, "plane": [1, 0, 0]}
        db_file (str): path to database file
        vasp_input_set (VaspInputSet): vasp_input_set corresponding to the grain boundary calculation
        parents (Fireworks or list of ints): parent FWs
        vasp_cmd (str): vasp command
        name (str): Name for the Firework

    Return:
        Firework
    """
    if gb_gen_params is None:
        raise ValueError("You are giving a bulk structure, please provide the "
                         "parameters to build grain boundary.")
    else:
        gb_gen_params.update(gb_gen_params)

    vis_for_transed_gb = vasp_input_set or MVLGBSet(bulk_structure, k_product=30)
    transformations = ["GrainBoundaryTransformation"]
    trans_params = [gb_gen_params]

    return TransmuterFW(structure=bulk_structure, transformations=transformations, name="gb_transmuter",
                        transformation_params=trans_params, copy_vasp_outputs=True, db_file=db_file,
                        vasp_cmd=vasp_cmd, parents=parents, vasp_input_set=vis_for_transed_gb)


def get_wf_gb_from_bulk(bulk_structure, gb_gen_params=None, tag=None, additional_info=None,
                        db_file=None, vasp_cmd="vasp"):
    """
    This is a workflow for grain boundary (gb) generation and calculation. If a bulk_structure and
     grain boundary (GB) generation parameters are specified, the workflow will relax the bulk structure
     first, then generate GB based on given parameters using TransmuterFW, which include the GB relaxation.



    Args:
        bulk_structure (Structure): bulk structure from which generate gbs after relaxation.
        gb_gen_params (dict): dictionary of gb generation parameters, used to generate gb. The
            details of description could be found in pymatgen.analysis.gb.gb import GBGenerator,
        tag (list): list of strings to tag the workflow, which will be inserted into the database
            make it easier to query data later. e.g.tag = ["mp-129"] to represent bcc-Mo wf.
        additional_info (Dict): the additional info of gb structure, which you want to add.
            The additional_info will be inserted into database.
        db_file (str): path to database file.
        vasp_cmd (str): vasp command

    Return:
        Workflow
    """
    fws, parents = [], []

    vis = MVLGBSet(bulk_structure, k_product=30)
    fws.append(OptimizeFW(structure=bulk_structure, vasp_input_set=vis,
                          vasp_cmd=vasp_cmd, db_file=db_file, name="bulk relaxation"))

    parents.append(fws[0])

    fws.append(get_gb_fw(bulk_structure=bulk_structure, gb_gen_params=gb_gen_params,
                         db_file=db_file, vasp_cmd=vasp_cmd, parents=parents))

    wf = Workflow(fws, name="{} gb workflow, e.g., {}".format(len(fws), fws[0].name))
    if additional_info:
        wf_additional_info = add_additional_fields_to_taskdocs(original_wf=wf, update_dict=additional_info)
        wf_with_tag_info = add_tags(wf_additional_info, tags_list=tag)
    else:
        wf_with_tag_info = add_tags(wf, tags_list=tag)

    return wf_with_tag_info


def get_wf_gb(gb, vasp_input_set=None, tag=None, additional_info=None, db_file=None, vasp_cmd="vasp"):
    """
     The workflow will directly run the relaxation for the given GB structure, while the related
     additional info could be added into the database as user specified.

     Args:
        gb (Structure/Gb): Grain boundary structure, could be Structure, or Gb (Gb object can
            be imported from pymatgen.analysis.gb.gb)
        vasp_input_set (VaspInputSet): vasp_input_set corresponding to the grain boundary calculation
        tag (list): list of strings to tag the workflow, which will be inserted into the database
            make it easier to query data later. e.g.tag = ["mp-129"] to represent bcc-Mo wf.:
        additional_info (Dict): the additional info of gb structure, which you want to add.
            The additional_info will be inserted into database.
        db_file (str): path to database file.
        vasp_cmd (str): vasp command

    Return:
        Workflow
    """
    fws, parents = [], []

    user_incar_settings = {"PREC": "Normal", "NPAR": 4, "ISMEAR": 1, "ENCUT": 400, "ICHARG": 2}
    vis_for_given_gb = vasp_input_set or MVLGBSet(gb, k_product=30,
                                                  user_incar_settings=user_incar_settings)

    name = gb.composition.reduced_formula
    fw = OptimizeFW(structure=gb, vasp_input_set=vis_for_given_gb, vasp_cmd=vasp_cmd,
                    parents=parents, db_file=db_file, name=name + "gb optimization")

    fws.append(fw)

    wf = Workflow(fws, name="{} gb workflow, e.g., {}".format(len(fws), fws[0].name))
    wf_additional_info = add_additional_fields_to_taskdocs(original_wf=wf, update_dict=additional_info)
    wf_with_tag_info = add_tags(wf_additional_info, tags_list=tag)
    return wf_with_tag_info
