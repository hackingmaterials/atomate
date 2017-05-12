# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines a workflow for adsorption on surfaces
"""

import numpy as np

from fireworks import Workflow

from atomate.utils.utils import get_logger
from atomate.vasp.fireworks.core import OptimizeFW, TransmuterFW, StaticFW

from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.core.surface import generate_all_slabs
from pymatgen.transformations.advanced_transformations import SlabTransformation
from pymatgen.transformations.standard_transformations import SupercellTransformation
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.sets import MVLSlabSet, MPRelaxSet
from pymatgen import Structure, Lattice

__author__ = 'Joseph Montoya'
__email__ = 'montoyjh@lbl.gov'

logger = get_logger(__name__)

# Default parameters for slabs (sip) and adsorbates (aip)
default_sip = {"ISIF": 0, "EDIFFG": -0.05}
default_aip = {"ISIF": 0, "AMIX": 0.1, "AMIX_MAG": 0.4, "BMIX": 0.0001,
               "BMIX_MAG": 0.0001, "POTIM": 0.3, "EDIFFG": -0.05}
default_slab_gen_params = {"max_index": 1, "min_slab_size": 7.0, "min_vacuum_size": 20.0,
                           "center_slab": True, "max_normal_search": None}


# TODO: @montoyjh - please fix the docstring -computron
# TODO: @montoyjh - this code is a huge mess. Please clean up, let me know when finished, and I'll review again. -computron

def get_wf_adsorption(structure, adsorbate_config, vasp_input_set=None, slab_gen_params=None,
                      vasp_cmd="vasp", db_file=None, conventional=True, slab_incar_params=None,
                      ads_incar_params=None, auto_dipole=True, use_bulk_coordination=False,
                      optimize_slab=True, optimize_bulk=True, find_ads_params=None):
    #TODO: add more details to docstring
    """
    Returns a workflow to calculate adsorption structures and surfaces.

    Firework 1 : write vasp input set for structural relaxation,
                 run vasp,
                 pass run location,
                 database insertion.

    Firework 2 - N: Optimize Slab and Adsorbate Structures
    
    Args:
        structure (Structure): input structure to be optimized and run
        adsorbate_config (dict): configuration dictionary, keys are strings
            corresponding to the miller index, values are molecules to be
            placed as adsorbates via the adsorption_site_finder
            e. g. {"111":Molecule("CO", [[0, 0, 0], [0, 0, 1.23]])}
        vasp_input_set (DictVaspInputSet): vasp input set.
        slab_gen_params
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.
        conventional
        slab_incar_params
        ads_incar_params
        auto_dipole
        use_bulk_coordination
        optimize_slab
        optimize_bulk
        find_ads_params

    Returns:
        Workflow
    """
    find_ads_params = find_ads_params or {}
    if conventional:
        sga = SpacegroupAnalyzer(structure)
        structure = sga.get_conventional_standard_structure()
    v = vasp_input_set or MVLSlabSet(structure, bulk=True)
    # Default VASP parameters for slabs and adsorbates,
    # may put this in a preset workflow eventually
    slab_incar_params = slab_incar_params or default_sip
    ads_incar_params = ads_incar_params or default_aip
    slab_gen_params = slab_gen_params or default_slab_gen_params
    if not v.incar.get("LDAU", None):
        ads_incar_params.update({"LDAU": False})
        slab_incar_params.update({"LDAU": False})

    fws = []
    if optimize_bulk:
        fws.append(OptimizeFW(structure=structure, vasp_input_set=v, vasp_cmd=vasp_cmd,
                              db_file=db_file, name="{}-structure optimization".format(
                structure.composition.reduced_formula)))

    max_index = max([int(i) for i in ''.join(adsorbate_config.keys())])  # TODO: @montoyjh - is this used anywhere? -computron
    slabs = generate_all_slabs(structure, **slab_gen_params)
    mi_strings = [''.join([str(i) for i in slab.miller_index])
                  for slab in slabs]
    for key in adsorbate_config.keys():
        if key not in mi_strings:
            raise ValueError("Miller index not in generated slab list. "
                             "Unique slabs are {}".format(mi_strings))
    
    for slab in slabs:
        mi_string = ''.join([str(i) for i in slab.miller_index])
        if mi_string in adsorbate_config.keys():
            # Add the slab optimize firework
            if auto_dipole:
                weights = np.array([site.species_and_occu.weight for site in slab])
                dipole_center = np.sum(weights*np.transpose(slab.frac_coords), axis=1)
                dipole_center /= np.sum(weights)
                dipole_dict = {"LDIPOL": "True",
                               "IDIPOL": 3,
                               "DIPOL": dipole_center}
                slab_incar_params.update(dipole_dict)
                ads_incar_params.update(dipole_dict)

            slab_trans_params = {"miller_index":slab.miller_index, "shift":slab.shift}
            slab_trans_params.update(slab_gen_params)
            slab_trans_params.pop("max_index")
            slab_trans = SlabTransformation(**slab_trans_params)  # TODO: @montoyjh - is this used anywhere? -computron

            slab.add_site_property("velocities", [[0., 0., 0.]]*slab.num_sites)  #  prevents problems with poscar/contcar conversion
            asf = AdsorbateSiteFinder(slab, selective_dynamics=True)  # TODO: @montoyjh - is this used anywhere? it's seemingly overridden below -computron
            if optimize_slab:
                fw_name = "{}_{} slab optimization".format(slab.composition.reduced_formula,
                                                           mi_string)
                fws.append(TransmuterFW(name=fw_name, structure=structure, db_file=db_file,
                                        vasp_cmd=vasp_cmd,
                                        transformations=["SlabTransformation",
                                                         "AddSitePropertyTransformation"],
                                        transformation_params=[
                                            slab_trans_params,
                                            {"site_properties": asf.slab.site_properties}],
                                        copy_vasp_outputs=True, parents=fws[0],
                                        vasp_input_set=MVLSlabSet(slab, user_incar_settings=slab_incar_params)))

            if use_bulk_coordination:
                asf = AdsorbateSiteFinder.from_bulk_and_miller(structure, slab.miller_index,
                                                               selective_dynamics=True)
            else:
                asf = AdsorbateSiteFinder(slab, selective_dynamics=True)

            # Generate adsorbate configurations and add fws to workflow
            for molecule in adsorbate_config[mi_string]:
                structures = asf.generate_adsorption_structures(molecule, repeat=[2, 2, 1],
                                                                find_args=find_ads_params)
                for n, struct in enumerate(structures):
                    # Sort structure because InsertSites sorts the structure
                    struct = struct.get_sorted_structure()

                    # TODO: @montoyjh - here is one reason the code is a mess. You define
                    # add_fw_name here but it's not needed until 10 lines down. Plus, you interrupt
                    # the flow of an actual code block regarding the variable "struct" to add this
                    # needless interjection. I cleaned up some of this but there is much more to
                    # go. Define things at the place they are needed -computron
                    ads_fw_name = "{}-{}_{} adsorbate optimization {}".format(
                        molecule.composition.formula, structure.composition.reduced_formula,
                        mi_string, n)

                    struct.add_site_property("velocities", [[0., 0., 0.]]*struct.num_sites)  # prevents problems with poscar/contcar conversion
                    trans_ads = ["SlabTransformation", "SupercellTransformation",
                                 "InsertSitesTransformation", "AddSitePropertyTransformation"]
                    trans_supercell = SupercellTransformation.from_scaling_factors(
                        round(struct.lattice.a / slab.lattice.a), 
                        round(struct.lattice.b / slab.lattice.b))
                    ads_sites = [site for site in struct if 
                                 site.properties["surface_properties"] == "adsorbate"]
                    trans_ads_params = [slab_trans_params,
                                        {"scaling_matrix": trans_supercell.scaling_matrix},
                                        {"species": [site.species_string for site in ads_sites],
                                         "coords": [site.frac_coords.tolist() # convert for proper serialization
                                                   for site in ads_sites]},
                                        {"site_properties": struct.site_properties}]
                    vis_ads = MVLSlabSet(structure, user_incar_settings=ads_incar_params)
                    fws.append(TransmuterFW(name=ads_fw_name,
                                            structure=structure,
                                            transformations=trans_ads,
                                            transformation_params=trans_ads_params,
                                            copy_vasp_outputs=True,
                                            db_file=db_file,
                                            vasp_cmd=vasp_cmd,
                                            parents=fws[0],
                                            vasp_input_set=vis_ads))
    wfname = "{}:{}".format(structure.composition.reduced_formula, "Adsorbate calculations")
    return Workflow(fws, name=wfname)


def get_wf_adsorption_from_slab(slab, molecules, vasp_input_set=None, vasp_cmd="vasp", db_file=None, 
                                auto_dipole=True, slab_incar_params=None, ads_incar_params=None, 
                                optimize_slab=False, find_ads_params=None, name=""):
    """
    This workflow function generates a workflow starting from a slab, rather than a bulk structure.
    It's intended use is for workflows that begin from optimized slabs or that begin from previous
    calculations that have already optimized structures.

    Args:
        slab (Slab or structure):
        molecules (Molecule):
        vasp_input_set
        vasp_cmd
        db_file
        auto_dipole
        slab_incar_params
        ads_incar_params
        optimize_slab
        find_ads_params
        name
    """
    find_ads_params = find_ads_params or {}
    slab_incar_params = slab_incar_params or default_sip
    ads_incar_params = ads_incar_params or default_aip
    if not name:
        name = slab.composition.reduced_formula
        if hasattr(slab, 'miller_index') and not name:
            name += "_"+"".join([str(i) for i in slab.miller_index])
    # TODO: @montoyjh - looks to repeat code from the previous workflow (get_wf_adsorption). Not
    # sure there is much "unique" code in this one. You could for example have the previous
    # workflow call this workflow and append the workflow from this one to the previous one. Then
    # get_wf_adsorption should also be much shorter. i.e., for each slab get_wf_adsorption, get the
    # wf from this class and append that WFlow. -computron
    if auto_dipole:
        weights = np.array([site.species_and_occu.weight for site in slab])
        dipole_center = np.sum(weights*np.transpose(slab.frac_coords), axis=1)
        dipole_center /= np.sum(weights)
        dipole_dict = {"LDIPOL": "True",
                       "IDIPOL": 3,
                       "DIPOL": dipole_center}
        slab_incar_params.update(dipole_dict)
        ads_incar_params.update(dipole_dict)

    fws = []
    if optimize_slab:
        slab_fw_name = "{} slab optimization".format(name)
        fws.append(OptimizeFW(structure=slab, vasp_input_set=vasp_input_set, vasp_cmd=vasp_cmd,
                              db_file=db_file, name=slab_fw_name))

    for molecule in molecules:
        adsorption_structures = AdsorbateSiteFinder(slab).generate_adsorption_structures(
            molecule, find_args=find_ads_params)
        for n, struct in enumerate(adsorption_structures):
            ads_fw_name = "{}-{} adsorbate optimization {}".format(
                molecule.composition.reduced_formula, name, n)
            fws.append(OptimizeFW(structure=slab, vasp_input_set=vasp_input_set, vasp_cmd=vasp_cmd,
                                  db_file=db_file, name=ads_fw_name))
    wfname = "{}:{}".format(name, "Adsorbate calculations")
    return Workflow(fws, name=wfname)

# TODO: @montoyjh - I imagine a common thing to do is to compute the molecule separately, the
# slab separately, then molecule + slab together. Is there a single Wflow to do all that? -computron

def get_wf_molecules(molecules, vasp_input_sets=None, vibrations=False, min_vacuum_size=15.0,
                     vasp_cmd="vasp", db_file=None):
    """
    Returns a workflow to calculate molecular energies as references for the
    surface workflow.

    Firework 1 : write vasp input set for structural relaxation,
                 run vasp,
                 pass run location,
                 database insertion.
    Args:
        molecules (list of molecules): input structure to be optimized and run
        vasp_input_sets (DictVaspInputSet): vasp input set.
        vibrations
        min_vacuum_size
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.

    Returns:
        Workflow
    """
    fws = []
    vasp_input_sets = vasp_input_sets or [None for _ in molecules]
    for molecule, vis in zip(molecules, vasp_input_sets):
        m_struct = Structure(Lattice.cubic(min_vacuum_size), molecule.species_and_occu,
                             molecule.cart_coords, coords_are_cartesian=True)
        m_struct.translate_sites(list(range(len(m_struct))),
                                 np.array([0.5]*3) - np.average(m_struct.frac_coords, axis=0))
        # concurrent with MVLSlabSet
        user_incar_settings = {"ENCUT": 400, "ISMEAR": 0, "IBRION": 2, "ISIF": 0,
                               "EDIFF": 1e-6, "EDIFFG": -0.01, "POTIM": 0.02}
        v = vis or MPRelaxSet(m_struct, user_incar_settings=user_incar_settings,
                              user_kpoints_settings={"grid_density": 1})
        # TODO: are you just trying to avoid double relaxation? You can still use OptimizeFW, just
        # change the job_type parameter... -computron
        # Use StaticFW to avoid double relaxation
        fws.append(StaticFW(structure=m_struct, vasp_input_set=v,
                            vasp_cmd=vasp_cmd, db_file=db_file))
        if vibrations:
            # Turn off symmetry because it screws up automatic k-points
            user_incar_settings.update({"IBRION": 5, "ISYM": 0})
            v = MPRelaxSet(m_struct, user_incar_settings=user_incar_settings,
                           user_kpoints_settings={"grid_density": 1})
            # TODO: why not just modify the static FW to allow for custom incar params? -computron
            # This is a bit of a hack.  Seems static fireworks don't cleanly allow
            # for custom incar parameters.
            fw = TransmuterFW(structure=m_struct, vasp_input_set=v,
                              copy_vasp_outputs=True, parents=fws[-1],
                              vasp_cmd=vasp_cmd, db_file=db_file)

            fws.append(fw)

    wfname = "{}".format("Molecule calculations")
    return Workflow(fws, name=wfname)
