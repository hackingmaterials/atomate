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


def get_wf_gibbs_free_energy(structure, vasp_input_set=None, vasp_cmd="vasp", scaling_matrices=None,
                             db_file=None, reciprocal_density=None):
    """
    Returns a workflow.
    """
    vis_relax = vasp_input_set or MPRelaxSet(structure, force_gamma=True)

    vis_static = MPStaticSet(structure, force_gamma=True, lepsilon=True)
    vis_static.incar["ISIF"] = 2
    vis_static.incar["ISTART"] = 1
    for key in ["MAGMOM", "LDAUU", "LDAUJ", "LDAUL"]:
        vis_static.incar.pop(key, None)

    if reciprocal_density:
        vis_relax.config_dict["KPOINTS"].update({"reciprocal_density": reciprocal_density})
        vis_relax = vis_relax.__class__.from_dict(vis_relax.as_dict())
        vis_static.reciprocal_density = reciprocal_density

    fws=[]

    # Structure optimization firework
    fws.append(OptimizeFW(structure=structure, vasp_input_set=vis_relax, vasp_cmd=vasp_cmd, db_file=db_file))

    if scaling_matrices:
        for scaling_matrix in scaling_matrices:
            fw = TransmuterFW(name="gibbs deformation", structure=structure,
                              transformations=["SupercellTransformation"],
                              transformation_params=[{"scaling_matrix": scaling_matrix}],
                              vasp_input_set=vis_static, copy_vasp_outputs=True, parents=fws[0],
                              vasp_cmd=vasp_cmd, db_file=db_file)
            fws.append(fw)

    wfname = "{}:{}".format(structure.composition.reduced_formula, "gibbs free energy")
    return Workflow(fws, name=wfname)
