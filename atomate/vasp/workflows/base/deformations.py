# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the deformation workflow: structure optimization followed by transmuter fireworks.
"""

from fireworks import Workflow

from atomate.utils.utils import get_logger
from atomate.vasp.firetasks.glue_tasks import PassStressStrainData
from atomate.vasp.fireworks.core import OptimizeFW, TransmuterFW

from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet
from pymatgen.analysis.elasticity import symmetry_reduce

__author__ = 'Kiran Mathew'
__credits__ = 'Joseph Montoya'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


# TODO: @kmathew - add missing docstring (e.g., optimize_structure, relax_deformed..)
def get_wf_deformations(structure, deformations, name="deformation", vasp_input_set=None,
                        lepsilon=False, vasp_cmd="vasp", db_file=None, user_kpoints_settings=None,
                        pass_stress_strain=False, tag="", relax_deformed=False, optimize_structure=True,
                        symmetry_reduction=False, metadata=None, pass_kpoints=False):
    """
    Returns a structure deformation workflow.

    Firework 1 : structural relaxation
    Firework 2 - len(deformations): Deform the optimized structure and run static calculations.


    Args:
        structure (Structure): input structure to be optimized and run
        deformations (list of 3x3 array-likes): list of deformations
        name (str): some appropriate name for the transmuter fireworks.
        vasp_input_set (DictVaspInputSet): vasp input set.
        lepsilon (bool): whether or not compute static dielectric constant/normal modes
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.
        user_kpoints_settings (dict): example: {"grid_density": 7000}
        pass_stress_strain (bool): if True, stress and strain will be parsed and passed on.
        tag (str): some unique string that will be appended to the names of the fireworks so that
            the data from those tagged fireworks can be queried later during the analysis.
        metadata (dict): meta data
        pass_kpoints (bool): if True, will ensure that deformations keep the same k-point mesh
            as the optimization firework

    Returns:
        Workflow
    """

    # TODO: @kmathew - why not use UUID to generate a unique tag for the user
    # if they don't feel like inventing one themselves? -computron
    fws, parents = [], []

    # TODO: @kmathew - I don't see the need for this option. Better if a user can just take an
    # OptimizeStructure workflow and chain it before this one? It's better if the workflows can
    # concentrate on what they are doing best, and encourage the user to chain things together. You
    # could create a preset workflow that chains it for the user, but I don't see a need to do it
    # here. Maybe I'm wrong though? -computron

    # input set for relaxation
    vis_relax = vasp_input_set or MPRelaxSet(structure, force_gamma=True)
    if optimize_structure:
        # Structure optimization firework
        fws = [OptimizeFW(structure=structure, vasp_input_set=vis_relax, vasp_cmd=vasp_cmd,
                          db_file=db_file, name="{} structure optimization".format(tag))]
        parents = fws[0]

    uis_static = vis_relax.user_incar_settings or {}
    uis_static = uis_static.copy()
    uis_static.update({"ISIF": 2, "ISTART":1})
    if relax_deformed:
        uis_static.update({"IBRION": 2, "NSW": 99})

    # ensure TransmuterFWs have same kpoints as undeformed structure if specified
    if pass_kpoints:
        uks_static = {"kpoints_object":vis_relax.kpoints} 
    else:
        uks_static = vis_relax.user_kpoints_settings
    
    # TODO: @kmathew - see my previous comment about chaining workflows -computron
    # static input set
    vis_static = MPStaticSet(structure, force_gamma=True, lepsilon=lepsilon,
                             user_kpoints_settings=uks_static,
                             user_incar_settings=uis_static)

    # Do symmetry reduction and get corresponding symmops if specified
    if symmetry_reduction:
        deformations = symmetry_reduce(deformations, structure)
        symmops = deformations.values()
    else:
        symmops = [None]*len(deformations)

    # Deformation fireworks with the task to extract and pass stress-strain appended to it.
    for n, deformation in enumerate(deformations):
        fw = TransmuterFW(name="{} {} {}".format(tag, name, n), structure=structure,
                          transformations=['DeformStructureTransformation'],
                          transformation_params=[{"deformation": deformation.tolist()}],
                          vasp_input_set=vis_static, copy_vasp_outputs=optimize_structure,
                          parents=parents, vasp_cmd=vasp_cmd, db_file=db_file)

        if pass_stress_strain:
            fw.tasks.append(PassStressStrainData(number=n, symmops=symmops[n],
                                                 deformation=deformation.tolist()))
        fws.append(fw)

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)

    return Workflow(fws, name=wfname, metadata=metadata)
