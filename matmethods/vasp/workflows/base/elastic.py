# coding: utf-8

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import numpy as np
from fireworks import Workflow
from fireworks.features.dupefinder import DupeFinderBase

from matmethods.vasp.fireworks.core import OptimizeFW, TransmuterFW
from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.analysis.elasticity.strain import IndependentStrain, Deformation
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pymatgen.transformations.standard_transformations import \
    DeformStructureTransformation

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

    Firework 2: Determine which deformations to run based on symmetry if needed
                TODO: How to turn off children selectively?

    Firework 3 - N: Optimize Deformed Structure

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

    # MP standards for elastic vasp inputs, kpoints
    v = vasp_input_set or MPRelaxSet(structure, force_gamma=True)
    """
                                     user_incar_settings = {"ENCUT": 700,
                                                            "EDIFF": 1e-6})
    v.config_dict['KPOINTS'].update({"grid_density": 7000})
    """
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

if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest

    structure=PymatgenTest.get_structure("Si")
    import pdb, traceback, sys
    try:
        wf=get_wf_elastic_constant(structure, vasp_cmd = ">>vasp_cmd<<",
                                   db_file = ">>db_file<<")
        from fireworks import LaunchPad
        lpad = LaunchPad.auto_load()
        lpad.add_wf(wf)
    except:
        type, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)
