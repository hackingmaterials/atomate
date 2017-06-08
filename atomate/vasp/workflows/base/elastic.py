# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the elastic workflow
"""

from pymatgen.analysis.elasticity.strain import Deformation
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.sets import MPStaticSet

from fireworks import Firework, Workflow

from atomate.utils.utils import get_logger
from atomate.vasp.workflows.base.deformations import get_wf_deformations
from atomate.vasp.firetasks.parse_outputs import ElasticTensorToDb

__author__ = 'Shyam Dwaraknath, Joseph Montoya'
__email__ = 'shyamd@lbl.gov, montoyjh@lbl.gov'

logger = get_logger(__name__)


def get_wf_elastic_constant(structure, norm_deformations=None, shear_deformations=None,
                            additional_deformations=None, vasp_input_set=None, vasp_cmd="vasp",
                            db_file=None, user_kpoints_settings=None, add_analysis_task=True,
                            conventional=True, copy_vasp_outputs=True, tag="elastic"):
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
        conventional (bool): If True the input structure will be converted to the conventional
            standard structure.
        copy_vasp_outputs (bool): whether or not copy the outputs from the previous calc(usually
            structure optimization) before running the deformed structure calculations(transmuter
            fireworks).

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

    vis_static = vasp_input_set or MPStaticSet(structure, force_gamma=True, lepsilon=False,
                                               user_kpoints_settings=user_kpoints_settings)

    wf_elastic = get_wf_deformations(structure, deformations, vasp_cmd=vasp_cmd, db_file=db_file,
                                     pass_stress_strain=True, name="deformation", tag=tag,
                                     vasp_input_set=vis_static, copy_vasp_outputs=copy_vasp_outputs)

    if add_analysis_task:
        fw_analysis = Firework(ElasticTensorToDb(structure=structure, db_file=db_file),
                               name="Analyze Elastic Data", spec={"_allow_fizzled_parents": True})
        wf_elastic.append_wf(Workflow.from_Firework(fw_analysis), wf_elastic.leaf_fw_ids)

    wf_elastic.name = "{}:{}".format(structure.composition.reduced_formula, "elastic constants")

    return wf_elastic


if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest

    structure = PymatgenTest.get_structure("Si")
    wf = get_wf_elastic_constant(structure)
