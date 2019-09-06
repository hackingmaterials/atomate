# coding: utf-8

"""
This module defines the elastic workflow
"""

import numpy as np

from pymatgen.analysis.elasticity.strain import Deformation, Strain
from pymatgen.core.tensors import symmetry_reduce
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.sets import MPStaticSet

from fireworks import Firework, Workflow

from atomate.utils.utils import get_logger, get_fws_and_tasks
from atomate.vasp.workflows.base.deformations import get_wf_deformations
from atomate.vasp.firetasks.parse_outputs import ElasticTensorToDb
from atomate.vasp.firetasks.glue_tasks import pass_vasp_result

__author__ = 'Shyam Dwaraknath, Joseph Montoya'
__email__ = 'shyamd@lbl.gov, montoyjh@lbl.gov'

logger = get_logger(__name__)


def get_wf_elastic_constant(structure, strain_states=None, stencils=None,
                            db_file=None,
                            conventional=False, order=2, vasp_input_set=None,
                            analysis=True,
                            sym_reduce=False, tag='elastic',
                            copy_vasp_outputs=False, **kwargs):
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
        strain_states (list of Voigt-notation strains): list of ratios of nonzero elements
            of Voigt-notation strain, e. g. [(1, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0), etc.].
        stencils (list of floats, or list of list of floats): values of strain to multiply
            by for each strain state, i. e. stencil for the perturbation along the strain
            state direction, e. g. [-0.01, -0.005, 0.005, 0.01].  If a list of lists,
            stencils must correspond to each strain state provided.
        db_file (str): path to file containing the database credentials.
        conventional (bool): flag to convert input structure to conventional structure,
            defaults to False.
        order (int): order of the tensor expansion to be determined.  Defaults to 2 and
            currently supports up to 3.
        vasp_input_set (VaspInputSet): vasp input set to be used.  Defaults to static
            set with ionic relaxation parameters set.  Take care if replacing this,
            default ensures that ionic relaxation is done and that stress is calculated
            for each vasp run.
        analysis (bool): flag to indicate whether analysis task should be added
            and stresses and strains passed to that task
        sym_reduce (bool): Whether or not to apply symmetry reductions
        tag (str):
        copy_vasp_outputs (bool): whether or not to copy previous vasp outputs.
        kwargs (keyword arguments): additional kwargs to be passed to get_wf_deformations

    Returns:
        Workflow
    """
    # Convert to conventional if specified
    if conventional:
        structure = SpacegroupAnalyzer(
            structure).get_conventional_standard_structure()

    uis_elastic = {"IBRION": 2, "NSW": 99, "ISIF": 2, "ISTART": 1,
                   "PREC": "High"}
    vis = vasp_input_set or MPStaticSet(structure,
                                        user_incar_settings=uis_elastic)
    strains = []
    if strain_states is None:
        strain_states = get_default_strain_states(order)
    if stencils is None:
        stencils = [np.linspace(-0.01, 0.01, 5 + (order - 2) * 2)] * len(
            strain_states)
    if np.array(stencils).ndim == 1:
        stencils = [stencils] * len(strain_states)
    for state, stencil in zip(strain_states, stencils):
        strains.extend(
            [Strain.from_voigt(s * np.array(state)) for s in stencil])

    # Remove zero strains
    strains = [strain for strain in strains if not (abs(strain) < 1e-10).all()]
    vstrains = [strain.voigt for strain in strains]
    if np.linalg.matrix_rank(vstrains) < 6:
        # TODO: check for sufficiency of input for nth order
        raise ValueError("Strain list is insufficient to fit an elastic tensor")

    deformations = [s.get_deformation_matrix() for s in strains]

    if sym_reduce:
        # Note this casts deformations to a TensorMapping
        # with unique deformations as keys to symmops
        deformations = symmetry_reduce(deformations, structure)

    wf_elastic = get_wf_deformations(structure, deformations, tag=tag,
                                     db_file=db_file,
                                     vasp_input_set=vis,
                                     copy_vasp_outputs=copy_vasp_outputs,
                                     **kwargs)
    if analysis:
        defo_fws_and_tasks = get_fws_and_tasks(wf_elastic,
                                               fw_name_constraint="deformation",
                                               task_name_constraint="Transmuted")
        for idx_fw, idx_t in defo_fws_and_tasks:
            defo = \
            wf_elastic.fws[idx_fw].tasks[idx_t]['transformation_params'][0][
                'deformation']
            pass_dict = {
                'strain': Deformation(defo).green_lagrange_strain.tolist(),
                'stress': '>>output.ionic_steps.-1.stress',
                'deformation_matrix': defo}
            if sym_reduce:
                pass_dict.update({'symmops': deformations[defo]})

            mod_spec_key = "deformation_tasks->{}".format(idx_fw)
            pass_task = pass_vasp_result(pass_dict=pass_dict,
                                         mod_spec_key=mod_spec_key)
            wf_elastic.fws[idx_fw].tasks.append(pass_task)

        fw_analysis = Firework(
            ElasticTensorToDb(structure=structure, db_file=db_file,
                              order=order, fw_spec_field='tags'),
            name="Analyze Elastic Data", spec={"_allow_fizzled_parents": True})
        wf_elastic.append_wf(Workflow.from_Firework(fw_analysis),
                             wf_elastic.leaf_fw_ids)

    wf_elastic.name = "{}:{}".format(structure.composition.reduced_formula,
                                     "elastic constants")

    return wf_elastic


def get_default_strain_states(order=2):
    """
    Generates a list of "strain-states"
    """
    inds = [(i,) for i in range(6)]
    # Note that these strain states may not be minimal
    if order > 2:
        inds.extend(
            [(0, i) for i in range(1, 5)] + [(1, 2), (3, 4), (3, 5), (4, 5)])
        if order > 3:
            inds.extend([(0, 1, 2), (0, 1, 3), (0, 1, 4), (0, 1, 5), (0, 2, 3),
                         (0, 2, 4), (0, 2, 5), (1, 2, 3), (1, 2, 4), (1, 2, 5),
                         (2, 3, 4), (2, 3, 5), (2, 4, 5), (3, 4, 5)])
            if order > 4:
                raise ValueError(
                    "Standard deformations for tensors higher than rank 4 not yet determined")
    strain_states = np.zeros((len(inds), 6))
    for n, i in enumerate(inds):
        np.put(strain_states[n], i, 1)
    strain_states[:, 3:] *= 2
    return strain_states.tolist()
