# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines a workflow for adsorption on surfaces
"""

import numpy as np

from fireworks import Workflow

from atomate.vasp.fireworks.core import OptimizeFW, TransmuterFW, StaticFW

from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.core.surface import generate_all_slabs, Slab
from pymatgen.transformations.advanced_transformations import SlabTransformation
from pymatgen.transformations.standard_transformations import SupercellTransformation
from pymatgen.io.vasp.sets import MVLSlabSet, MPRelaxSet
from pymatgen import Structure, Lattice

__author__ = 'Joseph Montoya'
__email__ = 'montoyjh@lbl.gov'


def get_slab_fw(slab, bulk_structure=None, slab_gen_params={}, db_file=None, vasp_input_set=None,
                parents=None, vasp_cmd="vasp", name=""):
    """
    Function to generate a a slab firework.  Returns a TransmuterFW if bulk_structure is specified,
    constructing the necessary transformations from the slab and slab generator parameters,
    or an OptimizeFW if only a slab is specified.

    Args:
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

    Returns:
        Firework
    """
    vasp_input_set = vasp_input_set or MVLSlabSet(slab)

    # If a bulk_structure is specified, generate the set of transformations, else
    # just create an optimize FW with the slab
    if bulk_structure:
        if not isinstance(slab, Slab):
            raise ValueError("structure input to get_slab_fw requires slab to be a slab object!")
        slab_trans_params = {"miller_index": slab.miller_index, "shift":slab.shift}
        slab_trans_params.update(slab_gen_params)

        # Get supercell parameters
        trans_struct = SlabTransformation(**slab_trans_params)
        slab_from_bulk = trans_struct.apply_transformation(bulk_structure)
        supercell_trans = SupercellTransformation.from_scaling_factors(
                round(slab.lattice.a / slab_from_bulk.lattice.a),
                round(slab.lattice.b / slab_from_bulk.lattice.b))

        # Get adsorbates for InsertSitesTransformation
        if "adsorbate" in slab.site_properties.get("surface_properties", [None]):
            ads_sites = [site for site in slab if site.properties["surface_properties"]=="adsorbate"]
        else:
            ads_sites = []
        transformations = ["SlabTransformation", "SupercellTransformation",
                           "InsertSitesTransformation", "AddSitePropertyTransformation"]
        trans_params = [slab_trans_params, {"scaling_matrix":supercell_trans.scaling_matrix},
                        {"species": [site.species_string for site in ads_sites],
                         "coords": [site.frac_coords for site in ads_sites]},
                        {"site_properties": slab.site_properties}]
        return TransmuterFW(name=name, structure=bulk_structure, transformations=transformations,
                            transformation_params=trans_params, copy_vasp_outputs=True, db_file=db_file,
                            vasp_cmd=vasp_cmd, parents=parents, vasp_input_set=vasp_input_set)
    else:
        return OptimizeFW(name=name, structure=slab, vasp_input_set=vasp_input_set, vasp_cmd=vasp_cmd,
                          db_file=db_file, parents=parents, job_type="normal")


def get_wf_surface(slabs, molecules=[], bulk_structure=None, slab_gen_params=None, vasp_cmd="vasp",
                   db_file=None, ads_structures_params={}, add_molecules_in_box=False):
    """

    Args:
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
        add_molecules_in_box (boolean): flag to add calculation of molecule
            energies to the workflow
        db_file (string): path to database file
        vasp_cmd (string): vasp command

    Returns:
        Workflow
    """
    fws, parents = [], []

    if bulk_structure:
        vis = MVLSlabSet(bulk_structure, bulk=True)
        fws.append(OptimizeFW(bulk_structure, vasp_input_set=vis, vasp_cmd="vasp", db_file=db_file))
        parents = fws[0]

    for slab in slabs:
        name = slab.composition.reduced_formula
        if getattr(slab, "miller_index", None):
            name += "_{}".format(slab.miller_index)
        fws.append(get_slab_fw(slab, bulk_structure, slab_gen_params, db_file=db_file, 
            vasp_cmd=vasp_cmd, parents=parents, name=name+" slab optimization"))
        for molecule in molecules:
            ads_slabs = AdsorbateSiteFinder(slab).generate_adsorption_structures(
                molecule, **ads_structures_params)
            for n, ads_slab in enumerate(ads_slabs):
                ads_name = "{}-{} adsorbate optimization {}".format(molecule.composition.formula, name, n)
                fws.append(get_slab_fw(ads_slab, bulk_structure, slab_gen_params, db_file=db_file,
                                       vasp_cmd=vasp_cmd, parents=parents, name=ads_name))
    if add_molecules_in_box:
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


def get_wf_surface_all_slabs(bulk_structure, molecules, max_index=1, slab_gen_params=None, **kwargs):
    """
    Convenience constructor that allows a user to construct a workflow
    that finds all adsorption configurations (or slabs) for a given
    max miller index.

    Args:
        bulk_structure (Structure): bulk structure from which to construct slabs
        molecules (list of Molecules): adsorbates to place on surfaces
        max_index (int): max miller index
        slab_gen_params (dict): dictionary of kwargs for generate_all_slabs

    Returns:
        Workflow
    """
    sgp = slab_gen_params or {"min_slab_size": 7.0, "min_vacuum_size": 20.0}
    slabs = generate_all_slabs(bulk_structure, max_index=max_index, **sgp)
    return get_wf_surface(slabs, molecules, bulk_structure, sgp, **kwargs)
