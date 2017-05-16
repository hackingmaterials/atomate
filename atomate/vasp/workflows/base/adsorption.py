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

# TODO: @montoyjh - please fix the docstring -computron

def SlabFW(slab, bulk_structure=None, slab_gen_params={}, vasp_input_set=None, 
           parents=None, db_file=None, vasp_cmd="vasp", name=""):
    """
    Constructor for a slab firework

    slab (Slab or Structure): structure or slab corresponding 
        to the slab to be calculated
    bulk_structure (Structure): bulk structure corresponding to slab, if
        provided, slab firework is constructed as a TransmuterFW using
        the necessary transformations to get the slab from the bulk
    slab_gen_params (dict): dictionary of slab generation parameters
        used to generate the slab, necessary to get the slab
        that corresponds to the bulk structure
    vasp_input_set (VaspInputSet): vasp_input_set corresponding to
        the slab calculation
    parents (Fireworks or list of ints): parent FWs
    db_file (string): path to database file
    vasp_cmd (string): vasp command
    """

    vasp_input_set = vasp_input_set or MVLSlabSet(slab)

    # Get adsorbate info and generate name
    if getattr(slab, "miller_index", None):
        name = "-{}".format(slab.miller_index) + name
    if "adsorbate" in slab.site_properties.get("surface_properties", [None]):
        ads_sites = [site for site in slab if 
                     site.properties['surface_properties'] == "adsorbate"]
        other_sites = [site for site in slab if 
                       site.properties['surface_properties'] != "adsorbate"]
        ads_name = Structure.from_sites(ads_sites).composition.reduced_formula_abc
        slab_name = Structure.from_sites(other_sites).composition.pretty_formula
        name = "adsorption {}_{}".format(ads_name, slab_name) + name
    else:
        ads_sites = []
        name = "slab {}".format(slab.composition.pretty_formula) + name

    if bulk_structure:
        if not isinstance(slab, Slab):
            raise ValueError("structure input to SlabFW requires slab to be a slab object!")
        slab_trans_params = {"miller_index": slab.miller_index, "shift":slab.shift}
        slab_trans_params.update(slab_gen_params)
        # Get supercell parameters
        trans_struct = SlabTransformation(**slab_trans_params)
        slab_from_bulk = trans_struct.apply_transformation(bulk_structure)
        supercell_trans = SupercellTransformation.from_scaling_factors(
                round(slab.lattice.a / slab_from_bulk.lattice.a),
                round(slab.lattice.b / slab_from_bulk.lattice.b))
        # Get adsorbates
                transformations = ["SlabTransformation", "SupercellTransformation",
                           "InsertSitesTransformation", "AddSitePropertyTransformation"]
        trans_params = [slab_trans_params, {"scaling_matrix":supercell_trans.scaling_matrix},
                        {"species": [site.species_string for site in ads_sites],
                         "coords": [site.frac_coords for site in ads_sites]}
                        {"site_properties": slab.site_properties}]
        # Designate name based on whether adsorbates or not
        return TransmuterFW(name=name, structure=structure, transformations=transformations,
                transformation_params=trans_params, copy_vasp_outputs=True, db_file=db_file,
                vasp_cmd=vasp_cmd, parents=parents, vasp_input_set=vasp_input_set)
    else:
        return OptimizeFW(structure=slab, vasp_input_set=vasp_input_set, vasp_cmd=vasp_cmd,
                db_file=db_file, vasp_cmd=vasp_cmd, parents=parents, job_type="normal")


def get_wf_surface(slabs, molecules=[], bulk_structure=None, slab_gen_params=None,
                   vasp_cmd="vasp", db_file=None, auto_dipole=False, ads_structures_params={},
                   vasp_input_set=None, add_molecules=False):
    """
    slabs (list of Slabs or Structures): slabs to calculate
    molecules (list of Molecules): molecules to place as adsorbates
    bulk_structure (Structure): bulk structure from which generate slabs
        after reoptimization.  If supplied, workflow will begin with
        bulk structure optimization.
    slab_gen_params (dict): dictionary of slab generation parameters
        used to generate the slab, necessary to get the slab
        that corresponds to the bulk structure if in that mode
    ads_structures_params (dict): parameters to be supplied as
        kwargs to AdsorbateSiteFinder.generate_adsorption_structures
    add_molecules (boolean): flag to add calculation of molecule
        energies to the workflow
    db_file (string): path to database file
    vasp_cmd (string): vasp command

    """
    vis = vasp_input_set or MVLSlabSet
    fws, parents = [], []
    if structure:
        fws.append(OptimizeFW(bulk_structure, vasp_cmd="vasp", db_file=None))
    for slab in slabs:
        fws.append(SlabFW(slab, bulk_structure, slab_gen_params, db_file=db_file, 
                          vasp_cmd=vasp_cmd, parents=parents))
        for molecule in molecules:
            ads_slabs = AdsorbateSiteFinder(slab).generate_adsorption_structures(
                    molecule, **ads_structures_params)
            for ads_slab in ads_slabs:
                fws.append(SlabFW(ads_slab, bulk_structure, slab_gen_params, db_file=db_file, 
                                  vasp_cmd=vasp_cmd, parents=parents))
    if add_molecules:
        for molecule in molecules:
            # molecule in box
            m_struct = Structure(Lattice.cubic(10), molecule.species_and_occu,
                                 molecule.cart_coords, coords_are_cartesian=True)
            m_struct.translate_sites(list(range(len(m_struct))),
                                 np.array([0.5]*3) - np.average(m_struct.frac_coords, axis=0))
            vis = MVLSlabSet(m_struct)
            fws.append(OptimizeFW(molecule, job_type="normal", vasp_input_set=vis,
                                  db_file=db_file, vasp_cmd=vasp_cmd))
    # TODO: add analysis framework
    return Workflow(fws, name="")

if __name__ == "__main__":
    pass
