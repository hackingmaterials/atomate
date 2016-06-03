# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines workflows for bandstructure calculations.
"""

from fireworks import Workflow
import os
from pymatgen.io.vasp.sets import MPVaspInputSet

from matmethods.vasp.fireworks.core import OptimizeFW, StaticFW, HSEBSFW, LcalcpolFW, NonSCFFW
from pymatgen.io.vasp.sets import MPStaticSet

from monty.serialization import loadfn

def generate_ferroelectric_template(nimages=5,optimize=True,band_hse=True):
    """
    This function returns a ferroelectric workflow for a given number of interpolations for the polarization calculation
    step.

    Args:
        nimages: number of interpolations for polarization calculations.
        optimize: optimize polar and nonpolar structures
        band_hse: calculate HSE band gap for polar structure

    Returns:
        YAML workflow

    """

    """
    CURRENT PROBLEMS WITH THIS FUNCTION

    Several of these Fireworks use previous jobs for the calculation. By default the calculation used for copying over
    the appropriate output files is the parent -- however the parent is not the appropriate previous calculation in
    several of these cases. Therefore, I have added the param 'prev' to several of the jobs.

    Also, uncertain how multiple parents should be attributed (should I be using links instead of 'parents'?).

    """

    yaml = ""

    if optimize:
        yaml= """fireworks:
- fw: matmethods.vasp.fireworks.core.OptimizeFW
    params:
        structure: 'polar'
        """

    yaml += """
- fw: matmethods.vasp.fireworks.core.StaticFW
    params:
        parents: 0
        structure: 'polar'
- fw: matmethods.vasp.fireworks.core.NonSCFFW
    params:
        parents: 1
        mode: uniform
        structure: 'polar'
- fw: matmethods.vasp.fireworks.core.NonSCFFW
    params:
        parents: 1
        mode: line
        structure: 'polar'
        defuse_children:
            band_gap:
                lt: 0 """

    if optimize:
        yaml += """
- fw: matmethods.vasp.fireworks.core.OptimizeFW
    params:
        parents: [2,3]
        prev: None
        structure: 'nonpolar' """


    fireworks = 2
    if optimize:
        fireworks += 2

    current_lcalcpol_parent = '[2,3]'
    if optimize:
        current_lcalcpol_parent = '4'

    for i in range(nimages):
        if i==nimages-1 and optimize:
            yaml += """
- fw: matmethods.vasp.fireworks.core.StaticFW
    params:
        parents: {p}
        prev: {prev_calc}
        structure: {i}
            """
            yaml = yaml.format(p='2,3', i=i, prev_calc=4)

        elif (i != nimages-1 and i != 0) or (i==nimages-1 and not optimize):
            yaml += """
- fw: matmethods.vasp.fireworks.core.StaticFW
    params:
        parents: {p}
        prev: None
        structure: {i}"""

            p = '[2,3]'
            if optimize:
                p = '4'
            yaml = yaml.format(p=p,i=i)
            fireworks += 1
            current_lcalcpol_parent = fireworks

        yaml += """
- fw: matmethods.vasp.fireworks.core.LcalcpolFW
    params:
        parents: {p}
        structure: {i}"""
        yaml = yaml.format(p=current_lcalcpol_parent,i=i)
        if i==0:
            yaml +="""
        prev: {prev_calc}"""
            yaml = yaml.format(prev_calc=1)
        fireworks += 1

# Also may want to add fw (or perhaps task) that calculates effective polarization.
# Alternatively, we can leave this to post-processing the database.

    if band_hse:
        yaml += """
- fw: matmethods.vasp.fireworks.core.HSEBSFW
    params:
        parents: [2,3]
        structure: 'polar' """

        fireworks += 1

    return yaml

def get_wf_ferroelectric(structure_polar, structure_nonpolar, nimages=5):
    """
    Simple workflow for ferroelectric search calculations. It includes these calculations:
        StaticFW x nimages
        LcalcpolFW x nimages
        NSCF x 1
        HSEBSFW x 1

    dependency of form:
        polar SCF
            polar NSCF
                interpolated and nonpolar SCF
                    interpolated and polar and nonpolar polarization
                polar HSEBSFW


    This does not yet include cuts or structure optimization. To incorporate such tasks we need to be able to:
        + pass interpolated structures to relevant Fireworks after relaxed polar and nonpolar structures
        + call completed job data or pass it such that cuts can be included as a FireTask

    Args:
        structure_polar (Structure): input polar structure
        structure_nonpolar (Structure):
            input nonpolar structure in polar symmetry setting. This is important since structures will be interpolated
            between each other

    Returns:
        Workflow
    """

    wfname = "{}-{}".format("FerroelectricSearch",
        structure_polar.composition.reduced_formula)


    commom_params={vasp_cmd: ">>vasp_cmd<<", db_file: '>>db_file<<'}

    # Interpolate structures from polar to nonpolar
    structs = structure_polar.interpolate(structure_nonpolar,nimages,interpolate_lattices=True)

    scf_fw = []
    polarization_fw = []

    scf_fw.append(StaticFW(s,name="scf_0",vasp_input_set=MPStaticSet(s),**commom_params))

    nscf = NonSCFFW(s[0],scf_fw[0],**commom_params)

    polarization_fw.append(
        LcalcpolFW(s, calc_loc="scf_0", parents=nscf, **commom_params))

    hsebs = HSEBSFW(s[0],nscf,**commom_params)

    for i,s in enumerate(structs[1:]):
        scf_fw.append(StaticFW(s,name="scf_"+str(i+1),parents=nscf,vasp_input_set=MPStaticSet(s),**commom_params))
        polarization_fw.append(
            LcalcpolFW(s,calc_loc="scf_"+str(i+1),parents=scf_fw[i+1],**commom_params))

    fws = scf_fw + polarization_fw + [nscf,hsebs]

    return Workflow(fws, name=wfname)

#if __name__ == "__main__":
#    from pymatgen.util.testing import PymatgenTest

# TO DO: replace this with BaTiO3 polar and nonpolar structures
#    structure = PymatgenTest.get_structure("Si")
#    wf = get_wf_ferroelectric(structure)
