# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the gibbs free energy workflow
"""

from fireworks import Workflow

from matmethods.utils.utils import get_logger
from matmethods.vasp.fireworks.core import OptimizeFW, TransmuterFW


from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


def get_wf_gibbs_free_energy(structure, vasp_input_set=None, vasp_cmd="vasp", deformations=None,
                             db_file=None, reciprocal_density=None):
    """
    Returns a workflow.
    """
    vis_relax = vasp_input_set or MPRelaxSet(structure, force_gamma=True)
    vis_static = MPStaticSet(structure, force_gamma=True, lepsilon=True)

    if reciprocal_density:
        vis_static.reciprocal_density = reciprocal_density

    fws = [OptimizeFW(structure=structure, vasp_input_set=vis_relax, vasp_cmd=vasp_cmd, db_file=db_file)]

    if deformations:
        for deform in deformations:
            fws.append(TransmuterFW(name="gibbs supercell transformation", structure=structure,
                                    transformations=['DeformStructureTransformation'],
                                    transformation_params=[{"deformation": deform.tolist()}],
                                    vasp_input_set=vis_static, copy_vasp_outputs=True, parents=fws[0],
                                    vasp_cmd=vasp_cmd, db_file=db_file))

    wfname = "{}:{}".format(structure.composition.reduced_formula, "gibbs free energy")
    return Workflow(fws, name=wfname)
