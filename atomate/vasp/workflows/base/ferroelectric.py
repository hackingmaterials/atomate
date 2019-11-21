# coding: utf-8


"""
This module defines the ferroelectric workflow
"""

from fireworks import Firework, Workflow

from atomate.utils.utils import get_logger, get_a_unique_id

from pymatgen import Structure

from atomate.vasp.fireworks.core import OptimizeFW
from atomate.vasp.fireworks.polarization import LcalcpolFW
from atomate.vasp.fireworks.core import HSEBSFW
from atomate.vasp.firetasks.parse_outputs import PolarizationToDb
from atomate.vasp.powerups import add_tags

__author__ = 'Tess Smidt'
__email__ = 'tsmidt@berkeley.edu'

logger = get_logger(__name__)


def get_wf_ferroelectric(polar_structure, nonpolar_structure, vasp_cmd="vasp", db_file=None,
                         vasp_input_set_polar="MPStaticSet", vasp_input_set_nonpolar="MPStaticSet",
                         relax=False, vasp_relax_input_set_polar=None, vasp_relax_input_set_nonpolar=None,
                         nimages=9, hse=False, add_analysis_task=False, wfid=None,
                         tags=None):
    """
    Returns a workflow to calculate the spontaneous polarization of polar_structure using
    a nonpolar reference phase structure and linear interpolations between the polar and
    nonpolar structure.

    The nonpolar and polar structures must be in the same space group setting and atoms ordered
    such that a linear interpolation can be performed to create intermediate structures along
    the distortion.

    For example, to calculate the polarization of orthorhombic BaTiO3 (space group 38) using
    the cubic structure (space group 221) as the nonpolar reference phase, we must transform
    the cubic to the orthorhombic setting. This can be accomplished using Bilbao Crystallographic
    Server's Structure Relations tool. (http://www.cryst.ehu.es/cryst/rel.html)

    Args:
        polar_structure (Structure): polar structure of candidate ferroelectric
        nonpolar_structure (Structure): nonpolar reference structure in polar setting
        vasp_input_set_polar (DictVaspInputSet): VASP polar input set. Defaults to MPStaticSet.
        vasp_input_set_nonpolar (DictVaspInputSet): VASP nonpolar input set. Defaults to MPStaticSet.
        vasp_relax_input_set_polar (DictVaspInputSet): VASP polar input set. Defaults to MPRelaxSet.
        vasp_relax_input_set_nonpolar (DictVaspInputSet): VASP nonpolar input set. Defaults to MPRelaxSet.
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.
        nimages: Number of interpolations calculated from polar to nonpolar structures, including the nonpolar.
            For example, nimages = 9 will calculate 8 interpolated structures. 8 interpolations + nonpolar = 9.
        add_analysis_task: Analyze polarization and energy trends as part of workflow. Default False.
        wfid (string): Unique worfklow id starting with "wfid_". If None this is atomatically generated (recommended).
        tags (list of strings): Additional tags to add such as identifiers for structures.

    Returns:

    """
    wf = []

    if wfid is None:
        wfid = 'wfid_' + get_a_unique_id()
    if tags is None:
        tags = []

    if relax:
        polar_relax = OptimizeFW(structure=polar_structure, name="_polar_relaxation",
                                 vasp_cmd=vasp_cmd, db_file=db_file, vasp_input_set=vasp_relax_input_set_polar)
        nonpolar_relax = OptimizeFW(structure=nonpolar_structure, name="_nonpolar_relaxation",
                                    vasp_cmd=vasp_cmd, db_file=db_file, vasp_input_set=vasp_relax_input_set_nonpolar)
        wf.append(polar_relax)
        wf.append(nonpolar_relax)
        parents_polar = polar_relax
        parents_nonpolar = nonpolar_relax
    else:
        parents_polar = None
        parents_nonpolar = None

    # Run polarization calculation on polar structure.
    # Defuse workflow if polar structure is metallic.
    polar = LcalcpolFW(structure=polar_structure,
                       name="_polar_polarization",
                       static_name="_polar_static",
                       parents=parents_polar,
                       vasp_cmd=vasp_cmd, db_file=db_file,
                       vasp_input_set=vasp_input_set_polar)

    # Run polarization calculation on nonpolar structure.
    # Defuse workflow if nonpolar structure is metallic.
    nonpolar = LcalcpolFW(structure=nonpolar_structure,
                          name="_nonpolar_polarization",
                          static_name="_nonpolar_static",
                          parents=parents_nonpolar,
                          vasp_cmd=vasp_cmd, db_file=db_file,
                          vasp_input_set=vasp_input_set_nonpolar)

    # Interpolation polarization
    interpolation = []
    # Interpolations start from one increment after polar and end prior to nonpolar.
    # The Structure.interpolate method adds an additional image for the nonpolar endpoint.
    # Defuse children fireworks if metallic.
    for i in range(1, nimages):
        # nonpolar_structure is being used as a dummy structure.
        # The structure will be replaced by the interpolated structure generated by
        # StaticInterpolatedFW.
        # Defuse workflow if interpolated structure is metallic.
        interpolation.append(
            LcalcpolFW(structure=polar_structure,
                       name="_interpolation_{}_polarization".format(str(i)),
                       static_name="_interpolation_{}_static".format(str(i)),
                       vasp_cmd=vasp_cmd, db_file=db_file,
                       vasp_input_set=vasp_input_set_polar, interpolate=True,
                       start="_polar_static",
                       end="_nonpolar_static",
                       nimages=nimages, this_image=i, parents=[polar, nonpolar]))

    wf.append(polar)
    wf.append(nonpolar)
    wf += interpolation

    # Add FireTask that uses Polarization object to store spontaneous polarization information
    if add_analysis_task:
        fw_analysis = Firework(PolarizationToDb(db_file=db_file),
                               parents=interpolation, name="_polarization_post_processing")
        wf.append(fw_analysis)

    # Run HSE band gap calculation
    if hse:
        # Run HSE calculation at band gap for polar calculation if polar structure is not metallic
        hse = HSEBSFW(structure=polar_structure, parents=polar, name="_polar_hse_gap", vasp_cmd=vasp_cmd,
                      db_file=db_file, calc_loc="_polar_polarization")
        wf.append(hse)

    # Create Workflow task and add tags to workflow
    workflow = Workflow(wf)
    workflow = add_tags(workflow, [wfid] + tags)

    return workflow
