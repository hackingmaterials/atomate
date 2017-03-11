# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the elastic workflow
"""
import itertools
import numpy as np

from pymatgen.analysis.elasticity.strain import Deformation, Strain
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen import Structure

from fireworks import Firework, Workflow

from atomate.utils.utils import get_logger, append_fw_wf
from atomate.vasp.workflows.base.deformations import get_wf_deformations
from atomate.vasp.firetasks.parse_outputs import ElasticTensorToDbTask, ToecToDbTask

__author__ = 'Shyam Dwaraknath, Joseph Montoya'
__email__ = 'shyamd@lbl.gov, montoyjh@lbl.gov'

logger = get_logger(__name__)


def get_wf_elastic_constant(structure, vasp_input_set=None, vasp_cmd="vasp", norm_deformations=None,
                            shear_deformations=None, additional_deformations=None, db_file=None,
                            user_kpoints_settings=None, conventional=True, optimize_structure=True,
                            symmetry_reduction=False):
    """
    Returns a workflow to calculate elastic constants.

    Firework 1 : write vasp input set for structural relaxation,
                 run vasp,
                 pass run location,
                 database insertion.

    Firework 2 - number of total deformations: Static runs on the deformed structures
    
    last Firework : Analyze Stress/Strain data and fit the elastic tensor

    Args:
        structure (Structure): input structure to be optimized and run.
        norm_deformations (list): list of values to for normal deformations.
        shear_deformations (list): list of values to for shear deformations.
        additional_deformations (list of 3x3 array-likes): list of additional deformations.
        vasp_input_set (DictVaspInputSet): vasp input set.
        vasp_cmd (str): command to run.
        db_file (str): path to file containing the database credentials.
        user_kpoints_settings (dict): example: {"grid_density": 7000}
        add_analysis_task (bool): boolean indicating whether to add analysis

    Returns:
        Workflow
    """
    # Convert to conventional
    if conventional:
        structure = SpacegroupAnalyzer(structure).get_conventional_standard_structure()
    # Generate deformations
    deformations = []

    if norm_deformations is not None:
        deformations.extend([Deformation.from_index_amount(ind, amount)
                             for ind in [(0, 0), (1, 1), (2, 2)]
                             for amount in norm_deformations])
    if shear_deformations is not None:
        deformations.extend([Deformation.from_index_amount(ind, amount)
                             for ind in [(0, 1), (0, 2), (1, 2)]
                             for amount in shear_deformations])

    if additional_deformations:
        deformations.extend([Deformation(defo_mat) for defo_mat in additional_deformations])

    if not deformations:
        raise ValueError("deformations list empty")

    wf_elastic = get_wf_deformations(structure, deformations, vasp_input_set=vasp_input_set, lepsilon=False, 
                                     vasp_cmd=vasp_cmd, db_file=db_file, user_kpoints_settings=user_kpoints_settings,
                                     pass_stress_strain=True, name="deformation", relax_deformed=True, tag="elastic", 
                                     symmetry_reduction=symmetry_reduction, optimize_structure=optimize_structure)

    fw_analysis = Firework(ElasticTensorToDbTask(structure=structure, db_file=db_file),
                           name="Analyze Elastic Data", spec={"_allow_fizzled_parents": True})
    append_fw_wf(wf_elastic, fw_analysis)

    wf_elastic.name = "{}:{}".format(structure.composition.reduced_formula, "elastic constants")

    return wf_elastic


def get_wf_toec(structure, vasp_input_set=None, vasp_cmd="vasp", db_file=None, 
                max_strain=0.05, stencil_res=7, indices=None, user_kpoints_settings=None, 
                conventional=True, optimize_structure=True, add_analysis_task=True):
    """
    Returns a workflow to calculate third-order elastic constants.

    Args:
        structure (Structure): input structure to be optimized and run. 
        vasp_input_set (DictVaspInputSet): vasp input set.
        vasp_cmd (str): command to run.
        db_file (str): path to file containing the database credentials.
        user_kpoints_settings (dict): example: {"grid_density": 7000}
        max_strain (float): maximum strain
        stencil_res (int): resolution on stencil to calculate second derivatives
        indices (list): list of indices e. g. [(1), (2), (3, 4)] to use for 
            strain states in deformed structures
        conventional (bool): flag to indicate whether to convert input structure 
            to conventional standard structure
        optimize_structure (bool): flag to indicate whether input structure
            should be optimized

    Returns:
        Workflow
    """
    # Convert to conventional
    if conventional:
        structure = SpacegroupAnalyzer(structure).get_conventional_standard_structure()
    # Generate deformations
    default_ind = [(i) for i in range(6)] + [(0, i) for i in range(1, 5)] \
            + [(1,2), (3,4), (3,5), (4,5)]
    indices = indices or default_ind
    strain_states = np.zeros((len(indices), 6))
    for n, index in enumerate(indices):
        np.put(strain_states[n], index, 1)
    strain_states[:, 3:] *= 2
    stencil = np.linspace(-max_strain, max_strain, stencil_res)
    stencil = stencil[np.nonzero(stencil)]
    deformations = [Strain.from_voigt(v*ss).deformation_matrix
                    for v, ss in itertools.product(stencil, strain_states)]

    wf_toec = get_wf_deformations(structure, deformations, vasp_input_set=vasp_input_set,
                                  lepsilon=False, vasp_cmd=vasp_cmd, db_file=db_file,
                                  user_kpoints_settings=user_kpoints_settings,
                                  pass_stress_strain=True, name="deformation",
                                  relax_deformed=True, tag="elastic",
                                  optimize_structure=optimize_structure)

    if add_analysis_task:
        fw_analysis = Firework(ToecToDbTask(structure=structure, db_file=db_file),
                               name="Analyze Elastic Data for TOEC",
                               spec={"_allow_fizzled_parents": True})
        append_fw_wf(wf_toec, fw_analysis)

    wf_toec.name = "{}:{}".format(structure.composition.reduced_formula, "third-order elastic constants")

    return wf_toec

if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest

    structure = PymatgenTest.get_structure("Si")
    #wf = get_wf_elastic_constant(structure)
    try:
        wf = get_wf_toec(structure)
    except:
        import sys, pdb, traceback
        type, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)
