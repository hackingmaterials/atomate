# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the elastic workflow
"""

from fireworks import Firework, Workflow

from matmethods.utils.utils import get_logger
from matmethods.vasp.fireworks.core import OptimizeFW, TransmuterFW
from matmethods.vasp.firetasks.glue_tasks import PassStressStrainData
from matmethods.vasp.firetasks.parse_outputs import ElasticTensorToDbTask

from pymatgen.analysis.elasticity.strain import Deformation
from pymatgen.io.vasp.sets import MPRelaxSet, DictSet, MPStaticSet
from pymatgen import Structure

__author__ = 'Shyam Dwaraknath, Joseph Montoya'
__email__ = 'shyamd@lbl.gov, montoyjh@lbl.gov'

logger = get_logger(__name__)


def get_wf_elastic_constant(structure, vasp_input_set=None, vasp_cmd="vasp", norm_deformations=None,
                            shear_deformations=None, additional_deformations=None, db_file=None,
                            reciprocal_density=None, add_analysis_task=True, lepsilon=False):
    """
    Returns a workflow to calculate elastic constants.

    Firework 1 : write vasp input set for structural relaxation,
                 run vasp,
                 pass run location,
                 database insertion.

    Firework 2 - 25: Optimize Deformed Structure
    
    Firework 26: Analyze Stress/Strain data and fit the elastic tensor

    Args:
        structure (Structure): input structure to be optimized and run
        norm_deformations (list): list of values to for normal deformations
        shear_deformations (list): list of values to for shear deformations
        additional_deformations (list of 3x3 array-likes): list of additional
            deformations to include
        vasp_input_set (DictVaspInputSet): vasp input set.
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.
        reciprocal_density (int): k-points per reciprocal atom by volume
        add_analysis_task (bool): boolean indicating whether to add analysis

    Returns:
        Workflow
    """
    vis_relax = vasp_input_set or MPRelaxSet(structure, force_gamma=True)

    vis_static = MPStaticSet(structure, force_gamma=True, lepsilon=lepsilon)
    vis_static.incar["ISIF"] = 2
    vis_static.incar["ISTART"] = 1
    for key in ["MAGMOM", "LDAUU", "LDAUJ", "LDAUL"]:
        vis_static.incar.pop(key, None)

    if reciprocal_density:
        vis_relax.config_dict["KPOINTS"].update({"reciprocal_density": reciprocal_density})
        vis_relax = vis_relax.__class__.from_dict(vis_relax.as_dict())
        vis_static.reciprocal_density = reciprocal_density

    fws=[]

    # Structure optimization firework
    fws.append(OptimizeFW(structure=structure, vasp_input_set=vis_relax, vasp_cmd=vasp_cmd, db_file=db_file))

    # Generate deformations
    deformations = []

    if norm_deformations:
        deformations.extend([Deformation.from_index_amount(ind, amount)
                             for ind in [(0, 0), (1, 1), (2, 2)]
                             for amount in norm_deformations])
    if shear_deformations:
        deformations.extend([Deformation.from_index_amount(ind, amount)
                             for ind in [(0, 1), (0, 2), (1, 2)]
                             for amount in shear_deformations])

    if additional_deformations:
        deformations.extend([Deformation(defo_mat) for defo_mat in additional_deformations])

    # Deformation fireworks with the task to extract and pass stress-strain appended to it.
    for deformation in deformations:
        fw = TransmuterFW(name="elastic deformation", structure=structure,
                          transformations=['DeformStructureTransformation'],
                          transformation_params=[{"deformation": deformation.tolist()}],
                          vasp_input_set=vis_static, copy_vasp_outputs=True, parents=fws[0],
                          vasp_cmd=vasp_cmd, db_file=db_file)
        fw.spec['_tasks'].append(PassStressStrainData(deformation=deformation.tolist()).to_dict())
        fws.append(fw)

    if add_analysis_task:
        fws.append(Firework(ElasticTensorToDbTask(structure=structure, db_file=db_file),
                            name="Analyze Elastic Data", parents=fws[1:],
                            spec={"_allow_fizzled_parents": True}))

    wfname = "{}:{}".format(structure.composition.reduced_formula, "elastic constants")
    return Workflow(fws, name=wfname)


if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest

    structure = PymatgenTest.get_structure("Si")
    wf = get_wf_elastic_constant(structure)
