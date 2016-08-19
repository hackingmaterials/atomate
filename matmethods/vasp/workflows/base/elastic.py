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
from pymatgen.io.vasp.sets import MPRelaxSet, DictSet
from pymatgen import Structure

__author__ = 'Shyam Dwaraknath, Joseph Montoya'
__email__ = 'shyamd@lbl.gov, montoyjh@lbl.gov'

logger = get_logger(__name__)


def get_wf_elastic_constant(structure, vasp_input_set=None, vasp_cmd="vasp",
                            norm_deformations=[-0.01, -0.005, 0.005, 0.01],
                            shear_deformations=[-0.06, -0.03, 0.03, 0.06],
                            db_file=None, reciprocal_density=None):
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
        vasp_input_set (DictVaspInputSet): vasp input set.
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.
        reciprocal_density (int): k-points per reciprocal atom by volume

    Returns:
        Workflow
    """

    v = vasp_input_set or MPRelaxSet(structure, force_gamma=True)
    if reciprocal_density:
        v.config_dict["KPOINTS"].update({"reciprocal_density": reciprocal_density})
        v = DictSet(structure, v.config_dict)
    fws=[]

    fws.append(OptimizeFW(structure=structure, vasp_input_set=v, vasp_cmd=vasp_cmd, db_file=db_file))

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

    def_incar_settings = v.incar.as_dict()
    def_incar_settings.update({"ISIF":2, "ISTART":1})
    for key in ["MAGMOM", "@module", "@class", "LDAUU", "LDAUJ", "LDAUL"]:
        def_incar_settings.pop(key, None)
    
    def_vasp_params = {"user_incar_settings":def_incar_settings}
    if reciprocal_density:
        def_vasp_params.update({"reciprocal_density":reciprocal_density})
    
    for deformation in deformations:
        # TODO: Maybe should be more general, needing to specify
        #   the vasp input set with a string is a bit unwieldy
        #   for complete customization of the INCAR parameters
        fw = TransmuterFW(name="elastic deformation",
                          structure=structure,
                          transformations=['DeformStructureTransformation'],
                          transformation_params=[{"deformation": deformation.tolist()}],
                          copy_vasp_outputs=True,
                          db_file=db_file,
                          vasp_cmd=vasp_cmd,
                          parents=fws[0],
                          vasp_input_params=def_vasp_params
                         )
        fw.spec['_tasks'].append(PassStressStrainData(deformation=deformation.tolist()).to_dict())
        fws.append(fw)
    
    fws.append(Firework(ElasticTensorToDbTask(structure=structure, db_file=db_file),
                        name="Analyze Elastic Data", parents=fws[1:],
                        spec={"_allow_fizzled_parents": True}))

    wfname = "{}:{}".format(structure.composition.reduced_formula, "elastic constants")
    return Workflow(fws, name=wfname)


if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest

    structure = PymatgenTest.get_structure("Si")
    wf = get_wf_elastic_constant(structure)
