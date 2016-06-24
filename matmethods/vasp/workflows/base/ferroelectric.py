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

        ideas for doing these things:

        re cuts:
            see MakeBandgapCut in parse_outputs.py
                This is a firetask that essentially redoes the work of VaspToDBTask but does not submit to the database.
                It uses the bandgap that gets inserted into the vasp database to make a cut.


        re passing interpolated structures:
            see GetInterpolatedPOSCAR vasp/firetasks/glue_tasks.py


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
    #structs = structure_polar.interpolate(structure_nonpolar,nimages,interpolate_lattices=True)

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
