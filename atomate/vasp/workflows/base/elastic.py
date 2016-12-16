# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the elastic workflow
"""

from pymatgen.analysis.elasticity.strain import Deformation
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen import Structure

from fireworks import Firework, Workflow

from atomate.utils.utils import get_logger, append_fw_wf
from atomate.vasp.workflows.base.deformations import get_wf_deformations
from atomate.vasp.firetasks.parse_outputs import ElasticTensorToDbTask

__author__ = 'Shyam Dwaraknath, Joseph Montoya'
__email__ = 'shyamd@lbl.gov, montoyjh@lbl.gov'

logger = get_logger(__name__)


def get_wf_elastic_constant(structure, vasp_input_set=None, vasp_cmd="vasp", norm_deformations=None,
                            shear_deformations=None, additional_deformations=None, db_file=None,
                            user_kpoints_settings=None, add_analysis_task=True, conventional=True):
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

    wf_elastic = get_wf_deformations(structure, deformations, vasp_input_set=vasp_input_set,
                                     lepsilon=False, vasp_cmd=vasp_cmd, db_file=db_file,
                                     user_kpoints_settings=user_kpoints_settings,
                                     pass_stress_strain=True, name="deformation",
                                     relax_deformed=True, tag="elastic")

    if add_analysis_task:
        fw_analysis = Firework(ElasticTensorToDbTask(structure=structure, db_file=db_file),
                               name="Analyze Elastic Data", spec={"_allow_fizzled_parents": True})
        append_fw_wf(wf_elastic, fw_analysis)

    wf_elastic.name = "{}:{}".format(structure.composition.reduced_formula, "elastic constants")

    return wf_elastic


if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest

    structure = PymatgenTest.get_structure("Si")
    wf = get_wf_elastic_constant(structure)
