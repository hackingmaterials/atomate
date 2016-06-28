# coding: utf-8

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from fireworks import Workflow

from matmethods.vasp.fireworks.core import OptimizeFW, TransmuterFW
from pymatgen.analysis.elasticity.strain import Deformation
from pymatgen.io.vasp.sets import MPRelaxSet

"""
This module defines the elastic workflow
"""

__author__ = 'Shyam Dwaraknath, Joseph Montoya'
__email__ = 'shyamd@lbl.gov, montoyjh@lbl.gov'

def get_wf_elastic_constant(structure, vasp_input_set=None, 
                            norm_deformations=[-0.01, -0.005, 0.005, 0.01],
                            shear_deformations=[-0.08, -0.04, 0.04, 0.08],
                            vasp_cmd="vasp", db_file=None):
    """
    Returns a workflow to calculate elastic constants.

    Firework 1 : write vasp input set for structural relaxation,
                 run vasp,
                 pass run location,
                 database insertion.

    Firework 2 - 25: Optimize Deformed Structure

    Args:
        structure (Structure): input structure to be optimized and run
        norm_deformations (list): list of values to for normal deformations
        shear_deformations (list): list of values to for shear deformations
        vasp_input_set (DictVaspInputSet): vasp input set.
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.

    Returns:
        Workflow
    """

    v = vasp_input_set or MPRelaxSet(structure, force_gamma=True)
    
    fws = []

    fws.append(OptimizeFW(structure=structure,
                          vasp_input_set=vasp_input_set,
                          vasp_cmd=vasp_cmd,
                          db_file=db_file))

    deformations = []
    # Generate deformations
    for ind in [(0, 0), (1, 1), (2, 2)]:
        for amount in norm_deformations:
            defo = Deformation.from_index_amount(ind, amount)
            deformations.append(defo)

    for ind in [(0, 1), (0, 2), (1, 2)]:
        for amount in shear_deformations:
            defo = Deformation.from_index_amount(ind, amount)
            deformations.append(defo)
    
    for deformation in deformations:
        fw = TransmuterFW(name="elastic deformation",
                          structure=structure,
                          vasp_input_set='MPRelaxSet',
                          transformations=['DeformStructureTransformation'],
                          transformation_params=[
                              {"deformation": deformation.tolist()}],
                          copy_vasp_outputs=True,
                          db_file=db_file,
                          vasp_cmd=vasp_cmd,
                          parents=fws[0],
                          vasp_input_params = {
                              "user_incar_settings":{"ISIF":2}}
                         )
        fws.append(fw)
    
    return Workflow(fws)
