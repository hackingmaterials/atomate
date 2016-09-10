# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the gibbs free energy workflow.
"""

from datetime import datetime

from fireworks import Firework, Workflow

from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet

from matmethods.utils.utils import get_logger
from matmethods.vasp.fireworks.core import OptimizeFW, TransmuterFW
from matmethods.vasp.firetasks.parse_outputs import GibbsFreeEnergyTask

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


def get_wf_gibbs_free_energy(structure, vasp_input_set=None, vasp_cmd="vasp", deformations=None,
                             db_file=None, reciprocal_density=None, t_step=10, t_min=0, t_max=1000,
                             mesh=(20, 20, 20), eos="vinet"):
    """
    Returns quasi-harmonic gibbs free energy workflow.
    Note: phonopy package is required for the final analysis step.

    Args:
        structure (Structure): input structure.
        vasp_input_set (VaspInputSet)
        vasp_cmd (str): vasp command to run.
        deformations (list): list of deformation matrices(list of lists).
        db_file (str): path to the db file.
        reciprocal_density (int)
        t_step (float): temperature step (in K)
        t_min (float): min temperature (in K)
        t_max (float): max temperature (in K)
        mesh (list/tuple): reciprocal space density
        eos (str): equation of state used for fitting the energies and the volumes.
            supported equation of states: vinet, murnaghan, birch_murnaghan

    Returns:
        Workflow
    """
    try:
        from phonopy import Phonopy
    except ImportError:
        logger.warn("'phonopy' package NOT installed. It is required for the final analysis step.")

    tag = datetime.utcnow().strftime('%Y-%m-%d-%H-%M-%S-%f')

    vis_relax = vasp_input_set or MPRelaxSet(structure, force_gamma=True)
    vis_static = MPStaticSet(structure, force_gamma=True, lepsilon=True)

    if reciprocal_density:
        vis_static.reciprocal_density = reciprocal_density

    fws = [OptimizeFW(structure=structure, vasp_input_set=vis_relax, vasp_cmd=vasp_cmd,
                      db_file=db_file, name="{} structure optimization".format(tag))]

    if deformations:
        for deform in deformations:
            fws.append(TransmuterFW(name="{} gibbs deformation".format(tag), structure=structure,
                                    transformations=['DeformStructureTransformation'],
                                    transformation_params=[{"deformation": deform}],
                                    vasp_input_set=vis_static, copy_vasp_outputs=True, parents=fws[0],
                                    vasp_cmd=vasp_cmd, db_file=db_file))

    fws.append(Firework(
        GibbsFreeEnergyTask(tag=tag, db_file=db_file, t_step=t_step, t_min=t_min, t_max=t_max,
                            mesh=mesh, eos=eos), name="Gibbs Free Energy", parents=fws[1:]))

    wfname = "{}:{}".format(structure.composition.reduced_formula, "gibbs free energy")
    return Workflow(fws, name=wfname)
