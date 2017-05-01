# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the elastic workflow
"""
import itertools
import numpy as np

from pymatgen.analysis.elasticity.strain import Deformation, Strain
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from fireworks import Firework, Workflow

from atomate.utils.utils import get_logger, append_fw_wf
from atomate.vasp.workflows.base.deformations import get_wf_deformations
from atomate.vasp.firetasks.parse_outputs import ElasticTensorToDbTask

# HACK
from atomate.vasp.powerups import add_common_powerups, add_modify_incar

__author__ = 'Shyam Dwaraknath, Joseph Montoya'
__email__ = 'shyamd@lbl.gov, montoyjh@lbl.gov'

logger = get_logger(__name__)


def get_wf_elastic_constant(structure, strain_states=None, stencils=None,
                            explicit_strains=[], db_file=None, conventional=True, 
                            order=2, **kwargs):
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

    strains = [Strain.from_voigt(strain) for strain in explicit_strains]
    if strain_states is None:
        strain_states = get_default_strain_states(order)
    if stencils is None:
        stencils = [np.linspace(-0.01, 0.01, 5 + (order-2)*2)]*len(strain_states)
    if np.array(stencils).ndim == 1:
        stencils = [stencils]*len(strain_states)
    for state, stencil in zip(strain_states, stencils):
        strains.extend([Strain.from_voigt(s*np.array(state)) for s in stencil])

    # Remove zero strains
    strains = [strain for strain in strains if not (abs(strain) < 1e-10).all()]
    vstrains = [strain.voigt for strain in strains]
    if np.linalg.matrix_rank(vstrains) < 6:
        # TODO: check for sufficiency of input for nth order
        raise ValueError("Strain list is insufficient to fit an elastic tensor")

    deformations = [s.deformation_matrix for s in strains]
    wf_elastic = get_wf_deformations(structure, deformations, pass_stress_strain=True, 
            name="deformation", relax_deformed=True, tag="elastic", db_file=db_file, **kwargs)

    fw_analysis = Firework(ElasticTensorToDbTask(structure=structure, db_file=db_file, order=order),
                           name="Analyze Elastic Data", spec={"_allow_fizzled_parents": True})
    append_fw_wf(wf_elastic, fw_analysis)

    wf_elastic.name = "{}:{}".format(structure.composition.reduced_formula, "elastic constants")

    return wf_elastic

def get_default_strain_states(order=2):
    """
    Generates a list of "strain-states"
    """
    inds = [(i,) for i in range(6)]
    if order > 2:
        inds.extend([(0, i) for i in range(1, 5)] + [(1,2), (3,4), (3,5), (4,5)])
        if order > 3:
            raise ValueError("Standard deformations for tensors higher than rank 4 not yet determined")
    strain_states = np.zeros((len(inds), 6))
    for n, i in enumerate(inds):
        np.put(strain_states[n], i, 1)
    strain_states[:, 3:] *= 2
    return strain_states.tolist()


def get_legacy_elastic_wf(structure, vasp_input_set=None, vasp_cmd=">>vasp_cmd<<", norm_deformations=None,
                            shear_deformations=None, additional_deformations=None, db_file=">>db_file<<",
                            user_kpoints_settings=None, add_analysis_task=True, conventional=True,
                            optimize_structure=True, symmetry_reduction=True, c=None):
    # Convert to conventional
    if conventional:
        structure = SpacegroupAnalyzer(structure).get_conventional_standard_structure()
    # Generate deformations
    user_kpoints_settings = user_kpoints_settings or {"grid_density": 7000}
    norm_deformations = norm_deformations or [0.005] #[-0.01, -0.005, 0.005, 0.01]
    shear_deformations = shear_deformations or [0.03] #[-0.06, -0.03, 0.03, 0.06]
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

    wf_elastic = get_wf_deformations(structure, deformations, vasp_input_set=vasp_input_set,
                                     lepsilon=False, vasp_cmd=vasp_cmd, db_file=db_file,
                                     user_kpoints_settings=user_kpoints_settings,
                                     pass_stress_strain=True, name="deformation",
                                     relax_deformed=True, tag="elastic",
                                     optimize_structure=optimize_structure, 
                                     symmetry_reduction=symmetry_reduction)

    if add_analysis_task:
        fw_analysis = Firework(ElasticTensorToDbTask(structure=structure, db_file=db_file),
                               name="Analyze Elastic Data", spec={"_allow_fizzled_parents": True})
        append_fw_wf(wf_elastic, fw_analysis)

    wf_elastic.name = "{}:{}".format(structure.composition.reduced_formula, "elastic constants")
    wf_elastic = add_modify_incar(wf_elastic, modify_incar_params = {"incar_update":{"ENCUT":700, "EDIFF":1e-6, "LAECHG":False}})
    wf_elastic = add_common_powerups(wf_elastic, c)

    return wf_elastic


if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest

    structure = PymatgenTest.get_structure("Si")
    #wf = get_wf_elastic_constant(structure)
    try:
        wf = get_wf_elastic_constant(structure, symmetry_reduction=True)
        wf2 = get_wf_elastic_constant(structure, order=3, symmetry_reduction=False)
    except:
        import sys, pdb, traceback
        type, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)
