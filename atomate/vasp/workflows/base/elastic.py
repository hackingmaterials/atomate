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
from atomate.vasp.firetasks.parse_outputs import ElasticTensorToDbTask

__author__ = 'Shyam Dwaraknath, Joseph Montoya'
__email__ = 'shyamd@lbl.gov, montoyjh@lbl.gov'

logger = get_logger(__name__)


def get_wf_elastic_constant(structure, norm_strains=None, shear_strains=None, 
                            additional_strains=None, db_file=None, conventional=True, 
                            toec=False, **kwargs):
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
        norm_strains (list of floats): list of values to for normal deformations.
        shear_strains (list of floats): list of values to for shear deformations.
        additional_strains (list of 3x3 array-likes): list of additional deformations.
        db_file (str): path to file containing the database credentials.
        toec (bool): whether to include TOEC analysis
        kwargs (keyword arguments): additional kwargs to be passed to get_wf_deformations

    Returns:
        Workflow
    """
    # Convert to conventional if specified
    if conventional:
        structure = SpacegroupAnalyzer(structure).get_conventional_standard_structure()

    # Generate strains from normal and shear values provided
    vstrains = []
    if norm_strains is not None:
        for norm_strain, ind in itertools.product(norm_strains, range(3)):
            vstrain = np.zeros(6)
            vstrain[ind] = norm_strain
            vstrains.append(vstrain)
    if shear_strains is not None:
        for shear_strain, ind in itertools.product(shear_strains, range(3,6)):
            vstrain = np.zeros(6)
            vstrain[ind] = shear_strain
            vstrains.append(vstrain)

    if additional_strains:
        vstrains.extend([Strain(strain_mat).voigt for strain_mat in additional_strains])

    if not vstrains or np.linalg.matrix_rank(np.array(vstrains)) < 6:
        raise ValueError("Strain list is insufficient to fit an elastic tensor")

    strains = [Strain.from_voigt(vs) for vs in vstrains]
    deformations = [s.deformation_matrix for s in strains]
    wf_elastic = get_wf_deformations(structure, deformations, pass_stress_strain=True, 
            name="deformation", relax_deformed=True, tag="elastic", **kwargs)

    fw_analysis = Firework(ElasticTensorToDbTask(structure=structure, db_file=db_file, toec=toec),
                           name="Analyze Elastic Data", spec={"_allow_fizzled_parents": True})
    append_fw_wf(wf_elastic, fw_analysis)

    wf_elastic.name = "{}:{}".format(structure.composition.reduced_formula, "elastic constants")

    return wf_elastic


def get_wf_toec(structure, max_strain=0.05, stencil_res=7, indices=None, **kwargs):
    """
    Returns a workflow to calculate third-order elastic constants.

    Args:
        structure (Structure): input structure to be optimized and run. 
        max_strain (float): maximum strain
        stencil_res (int): resolution on stencil to calculate second derivatives
        indices (list): list of indices e. g. [(1), (2), (3, 4)] to use for 
            strain states in deformed structures
        **kwargs (keyword arguments): kwargs to be passed to get_wf_elastic
    Returns:
        Workflow
    """
    if stencil_res % 2 != 1 or stencil_res < 5:
        raise ValueError("Stencil resolution for TOECs must be an odd integer greater than 5")

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
    strains = [Strain.from_voigt(v*ss) for v, ss in itertools.product(stencil, strain_states)]

    wf_toec = get_wf_elastic_constant(structure, norm_strains=[], shear_strains=[],
                                      additional_strains=strains, toec=True, **kwargs)
 
    wf_toec.name = "{}:{}".format(structure.composition.reduced_formula, "third-order elastic constants")

    return wf_toec

if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest

    structure = PymatgenTest.get_structure("Si")
    #wf = get_wf_elastic_constant(structure)
    try:
        wf = get_wf_elastic_constant(structure, norm_deformations=[0.01], 
                shear_deformations=[0.03], symmetry_reduction=True)
    except:
        import sys, pdb, traceback
        type, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)
