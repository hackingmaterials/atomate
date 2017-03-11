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


def get_wf_elastic_constant(structure, norm_deformations=None, shear_deformations=None, 
                            additional_deformations=None, db_file=None, conventional=True, 
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
        norm_deformations (list): list of values to for normal deformations.
        shear_deformations (list): list of values to for shear deformations.
        additional_deformations (list of 3x3 array-likes): list of additional deformations.
        db_file (str): path to file containing the database credentials.
        toec (bool): whether to include TOEC analysis
        kwargs (keyword arguments): additional kwargs to be passed to get_wf_deformations

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
    deformations = [Strain.from_voigt(v*ss).deformation_matrix
                    for v, ss in itertools.product(stencil, strain_states)]

    wf_toec = get_wf_elastic_constant(structure, norm_deformations=[], shear_deformations=[],
                                      additional_deformations=deformations, toec=True, **kwargs)
 
    wf_toec.name = "{}:{}".format(structure.composition.reduced_formula, "third-order elastic constants")

    return wf_toec

if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest

    structure = PymatgenTest.get_structure("Si")
    #wf = get_wf_elastic_constant(structure)
    try:
        wf = get_wf_elastic_constant(structure, norm_deformations=[0.01], 
                shear_deformations=[0.03], symmetry_reduction=True)
        from monty.serialization import dumpfn
        dumpfn(wf.to_dict(), "wf.json")
    except:
        import sys, pdb, traceback
        type, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)
