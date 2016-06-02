# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines workflows for bandstructure calculations.
"""

from fireworks import Workflow
import os
from pymatgen.io.vasp.sets import MPVaspInputSet

from monty.serialization import loadfn
from matmethods.utils.loaders import get_wf_from_spec_dict


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
    the appropriate output files is the parent -- however the parent is not the approciate previous calculation in
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
        parents: 2,3
        prev: None
        structure: 'nonpolar' """


    fireworks = 2
    if optimize:
        fireworks += 2

    current_lcalcpol_parent = 1

    for i in range(nimages):
        if i==nimages-1 and optimize:
            yaml += """
- fw: matmethods.vasp.fireworks.core.StaticFW
    params:
        parents: {p}
        prev: {prev_calc}
        structure: {i}
            """
            yaml = yaml.format(p='2,3', i=i, prev_calc=fireworks)

        elif (i != nimages-1 and i != 0) or (i==0 and not optimize):
            yaml += """
- fw: matmethods.vasp.fireworks.core.StaticFW
    params:
        parents: {p}
        prev: None
        structure: {i}"""

            yaml = yaml.format(p='2,3',i=i)
            fireworks += 1
            current_lcalcpol_parent = 4

        yaml += """
- fw: matmethods.vasp.fireworks.core.LcalcpolFW
    params:
        parents: {p}
        structure: {i}"""
        yaml = yaml.format(p=current_lcalcpol_parent,i=i)
        fireworks += 1

# Also may want to add fw (or perhaps task) that calculates effective polarization.
# Alternatively, we can leave this to post-processing the database.

    if band_hse:
        yaml += """
- fw: matmethods.vasp.fireworks.core.HSEBSFW
    params:
        parents: 2,3
        structure: 'polar' """

        fireworks += 1

    return yaml

def get_wf_ferroelectric(structure_polar, structure_nonpolar, vasp_input_set=None, vasp_cmd="vasp",
                         db_file=None):
    """
    NEED TO REWRITE TO DESCRIBE YAML GENERATION

    Return vasp workflow consisting of 4 fireworks:

    Firework 1 : write vasp input set for structural relaxation,
                 run vasp,
                 pass run location,
                 database insertion.

    Firework 2 : copy files from previous run,
                 write vasp input set for static run,
                 run vasp,
                 pass run location
                 database insertion.

    Firework 3 : copy files from previous run,
                 write vasp input set for non self-consistent
                 (constant charge density) run in
                 uniform mode,
                 run vasp,
                 pass run location
                 database insertion.

    Firework 4 : copy files from previous run,
                 write vasp input set for non self-consistent
                 (constant charge density) run in
                 line mode,
                 run vasp,
                 pass run location
                 database insertion.

    Firework 5 : (optional) HSE gap run

    Args:
        structure_polar (Structure): input polar structure
        structure_nonpolar (Structure):
            input nonpolar structure in polar symmetry setting. This is important since structures will be interpolated
            between each other
        vasp_input_set (DictVaspInputSet): vasp input set.
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.
        add_hse_gap (bool): add an HSE gap

    Returns:
        Workflow
    """

#    wf_file = "ferroelectric.yaml"

#    wfspec = loadfn(os.path.join(os.path.dirname(__file__), wf_file))

    wfspec = generate_ferroelectric_template()

    v = vasp_input_set or MPVaspInputSet()
    # We need to add lanthanide +U to the standard MP Vasp Inputs and copy as the FerroelectricSearchInputSet
    wfspec["fireworks"][0]["params"] = {"vasp_input_set": v.as_dict()}

    wfspec["common_params"] = {
        "vasp_cmd": vasp_cmd,
        "db_file": db_file
    }

    # can't use the default method because we need to specify multiple structures.
    #return get_wf_from_spec_dict(structure, wfspec)

    #### COPIED FROM loaders.get_wf_from_spec_dict

    dec = MontyDecoder()

    def load_class(dotpath):
        modname, classname = dotpath.rsplit(".", 1)
        mod = __import__(modname, globals(), locals(), [classname], 0)
        return getattr(mod, classname)

    def process_params(d):
        decoded = {}
        for k, v in d.items():
            if k.startswith("$"):
                if isinstance(v, list):
                    v = [os.path.expandvars(i) for i in v]
                elif isinstance(v, dict):
                    v = {k2: os.path.expandvars(v2) for k2, v2 in v.items()}
                else:
                    v = os.path.expandvars(v)
            decoded[k.strip("$")] = dec.process_decoded(v)
        return decoded

    fws = []
    common_params = process_params(wfspec.get("common_params", {}))

    # Note that interpolation of structures has to happen after optimizations


    for d in wfspec["fireworks"]:
        cls_ = load_class(d["fw"])
        params = process_params(d.get("params", {}))
        params.update(common_params)
        if "parents" in params:
            if isinstance(params["parents"], int):
                params["parents"] = fws[params["parents"]]
            else:
                p = []
                for parent_idx in params["parents"]:
                    p.append(fws[parent_idx])
                params["parents"] = p
        fws.append(cls_(structure, **params))


    wfname = "{}:{}".format(structure.composition.reduced_formula, wfspec["name"]) if \
        wfspec.get("name") else structure.composition.reduced_formula
    return Workflow(fws, name=wfname)

#if __name__ == "__main__":
#    from pymatgen.util.testing import PymatgenTest

# TO DO: replace this with BaTiO3 polar and nonpolar structures
#    structure = PymatgenTest.get_structure("Si")
#    wf = get_wf_ferroelectric(structure)
