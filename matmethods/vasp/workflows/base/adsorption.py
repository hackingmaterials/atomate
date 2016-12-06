# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines a workflow for adsorption on surfaces
"""

import json
from decimal import Decimal

import numpy as np

from fireworks import FireTaskBase, Firework, FWAction, Workflow
from fireworks.utilities.fw_serializers import DATETIME_HANDLER
from fireworks.utilities.fw_utilities import explicit_serialize

from matmethods.utils.utils import env_chk, get_logger
from matmethods.vasp.database import MMDb
from matmethods.vasp.fireworks.core import OptimizeFW, TransmuterFW

from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.core.surface import generate_all_slabs
from pymatgen.transformations.advanced_transformations import SlabTransformation
from pymatgen.transformations.standard_transformations import SupercellTransformation
from pymatgen.transformations.site_transformations import InsertSitesTransformation
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.sets import MVLSlabSet, MPRelaxSet, DictSet
from pymatgen import Structure, Lattice

__author__ = 'Joseph Montoya'
__email__ = 'montoyjh@lbl.gov'

logger = get_logger(__name__)

def get_wf_adsorption(structure, adsorbate_config, vasp_input_set=None,
                      min_slab_size=7.0, min_vacuum_size=20.0, center_slab=True,
                      max_normal_search=None, vasp_cmd="vasp", db_file=None, 
                      conventional=True, slab_incar_params=None, 
                      ads_incar_params=None, auto_dipole=True,
                      use_bulk_coordination=False):
    #TODO: add more details to docstring
    #TODO: this could use a general refactoring in the context of new ASF
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
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.

    Returns:
        Workflow
    """
    if conventional:
        sga = SpacegroupAnalyzer(structure)
        structure = sga.get_conventional_standard_structure()
    v = vasp_input_set or MVLSlabSet(structure, bulk=True)
    # Default VASP parameters for slabs and adsorbates,
    # may want to put this in a preset workflow
    slab_incar_params = slab_incar_params or {"ISIF":0}
    ads_incar_params = ads_incar_params or {"ISIF":0,
                                            "AMIX":0.1,
                                            "AMIX_MAG":0.4,
                                            "BMIX":0.0001,
                                            "BMIX_MAG":0.0001,
                                            "POTIM":0.25,
                                            "EDIFFG":-0.05}

    if not v.incar.get("LDAU", None):
        ads_incar_params.update({"LDAU":False})
        slab_incar_params.update({"LDAU":False})
    fws = []
    fws.append(OptimizeFW(structure=structure, vasp_input_set=v,
                          vasp_cmd=vasp_cmd, db_file=db_file))

    max_index = max([int(i) for i in ''.join(adsorbate_config.keys())])
    slabs = generate_all_slabs(structure, max_index=max_index,
                               min_slab_size=min_slab_size,
                               min_vacuum_size=min_vacuum_size, 
                               max_normal_search=max_normal_search,
                               center_slab=center_slab)
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
                weights = [site.species_and_occu.weight for site in slab]
                dipol_center = [0, 0, 0.5]
                #np.average(slab.frac_coords, axis = 1).tolist()
                dipole_dict = {"LDIPOL":"True",
                               "IDIPOL": 3,
                               "DIPOL": dipol_center}
                slab_incar_params.update(dipole_dict)
                ads_incar_params.update(dipole_dict)
            slab_trans_params = {"miller_index":slab.miller_index,
                           "min_slab_size":min_slab_size,
                           "min_vacuum_size":min_vacuum_size,
                           "shift":slab.shift,
                           "center_slab":True,
                           "max_normal_search":max_normal_search}
            slab_trans = SlabTransformation(**slab_trans_params)
            # TODO: name these more intelligently
            fw_name = "{}_{} slab optimization".format(
                slab.composition.reduced_formula, mi_string)
            vis_slab = MVLSlabSet(slab, user_incar_settings=slab_incar_params)
            fws.append(TransmuterFW(name=fw_name,
                                    structure = structure,
                                    transformations = ["SlabTransformation"],
                                    transformation_params=[slab_trans_params],
                                    copy_vasp_outputs=True,
                                    db_file=db_file,
                                    vasp_cmd=vasp_cmd,
                                    parents=fws[0],
                                    vasp_input_set = vis_slab)
                      )
            # Generate adsorbate configurations and add fws to workflow
            if use_bulk_coordination:
                asf = AdsorbateSiteFinder.from_bulk_and_miller(structure, slab.miller_index,
                                                               selective_dynamics=True)
            else:
                asf = AdsorbateSiteFinder(slab, selective_dynamics=True)
            for molecule in adsorbate_config[mi_string]:
                structures = asf.generate_adsorption_structures(molecule, repeat=[2, 2, 1])
                for struct in structures:
                    # Sort structure because InsertSites sorts the structure
                    struct = struct.get_sorted_structure()
                    ads_fw_name = "{}-{}_{} adsorbate optimization".format(
                        molecule.composition.formula,
                        structure.composition.reduced_formula, mi_string)
                    # Add velocities to avoid problems with poscar/contcar conversion
                    struct.add_site_property("velocities", [[0., 0., 0.]]*struct.num_sites)
                    trans_ads = ["SlabTransformation", "SupercellTransformation", 
                                  "InsertSitesTransformation", "AddSitePropertyTransformation"]
                    trans_supercell = SupercellTransformation.from_scaling_factors(
                        round(struct.lattice.a / slab.lattice.a), 
                        round(struct.lattice.b / slab.lattice.b))
                    ads_sites = [site for site in struct if 
                                 site.properties["surface_properties"]=="adsorbate"]
                    trans_ads_params = [slab_trans_params,
                                        {"scaling_matrix":trans_supercell.scaling_matrix},
                                        {"species":[site.species_string for site in ads_sites],
                                         "coords":[site.frac_coords.tolist() # convert for proper serialization
                                                   for site in ads_sites]},
                                        {"site_properties": struct.site_properties}]
                    vis_ads = MVLSlabSet(structure, user_incar_settings = ads_incar_params)
                    fws.append(TransmuterFW(name=ads_fw_name,
                                            structure = structure,
                                            transformations = trans_ads,
                                            transformation_params = trans_ads_params,
                                            copy_vasp_outputs=True,
                                            db_file=db_file,
                                            vasp_cmd=vasp_cmd,
                                            parents=fws[0],
                                            vasp_input_set = vis_ads))
    wfname = "{}:{}".format(structure.composition.reduced_formula, "Adsorbate calculations")
    return Workflow(fws, name=wfname)

def get_wf_molecules(molecules, vasp_input_sets=None,
                     min_vacuum_size=15.0, vasp_cmd="vasp",
                     db_file=None):
    """
    Returns a workflow to calculate molecular energies as references for the
    surface workflow.

    Firework 1 : write vasp input set for structural relaxation,
                 run vasp,
                 pass run location,
                 database insertion.
    Args:
        molecules (list of molecules): input structure to be optimized and run
        vasp_input_set (DictVaspInputSet): vasp input set.
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.

    Returns:
        Workflow
    """
    fws = []
    vasp_input_sets = vasp_input_sets or [None for m in molecules]
    for molecule, vis in zip(molecules, vasp_input_sets):
        m_struct = Structure(Lattice.cubic(min_vacuum_size), molecule.species_and_occu,
                             molecule.cart_coords, coords_are_cartesian=True)
        m_struct.translate_sites(list(range(len(m_struct))),
                                 np.array([0.5]*3) - np.average(m_struct.frac_coords,axis=0))
        v = vis or MPRelaxSet(m_struct, user_incar_settings={"ISMEAR":0, 
                                                             "IBRION":2, 
                                                             "ISIF":0,
                                                             "EDIFF":1e-6,
                                                             "EDIFFG":-0.01,
                                                             "POTIM":0.1},
                              user_kpoints_settings = {"grid_density":1}) #TODO think about this
        # There should also be a method to calculate hessian
        fws.append(OptimizeFW(structure=m_struct, vasp_input_set=v,
                              vasp_cmd=vasp_cmd, db_file=db_file))
    wfname = "{}".format("Molecule calculations")
    return Workflow(fws, name=wfname)


if __name__ == "__main__":
    from fireworks import LaunchPad
    lpad = LaunchPad.auto_load()
    from pymatgen.util.testing import PymatgenTest
    from pymatgen import Molecule, MPRester
    mpr = MPRester()
    pd = mpr.get_structures("mp-2")[0]
    co = mpr.get_structures("mp-54")[0]
    fe = mpr.get_structures("mp-13")[0]
    h = Molecule("H", [[0, 0, 0]])
    adsorbate_config = {k:[h] for k in ["111", "100", "110"]}
    structure = PymatgenTest.get_structure("Si")
    wf1 = get_wf_adsorption(pd, adsorbate_config)
    wf2 = get_wf_adsorption(co, {"001":[h]})
    wf3 = get_wf_adsorption(fe, {"100":[h],
                                 "110":[h]})
    fw_names = [fw.name for fw in wf1.fws]
    fw_names2 = [fw.name for fw in wf2.fws]
    fw_names3 = [fw.name for fw in wf3.fws]
    import pdb; pdb.set_trace()
    #wf2 = get_wf_molecules([co])
