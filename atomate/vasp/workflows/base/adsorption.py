# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines a workflow for adsorption on surfaces
"""

import numpy as np

from fireworks import Workflow

from atomate.vasp.fireworks.core import OptimizeFW, TransmuterFW

from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.core.surface import generate_all_slabs, Slab
from pymatgen.transformations.advanced_transformations import SlabTransformation
from pymatgen.transformations.standard_transformations import SupercellTransformation
from pymatgen.io.vasp.sets import MVLSlabSet
from pymatgen import Structure

__author__ = 'Joseph Montoya'
__email__ = 'montoyjh@lbl.gov'


#TODO: This doesn't currently allow for reconstructions, AFAIK
def get_slab_fw(slab, transmuter=False, slab_gen_params=None, db_file=None,
                vasp_input_set=None, parents=None, vasp_cmd="vasp", name=""):
    """
    Function to generate a a slab firework.  Returns a TransmuterFW if
    bulk_structure is specified, constructing the necessary transformations
    from the slab and slab generator parameters, or an OptimizeFW if only a
    slab is specified.

    Args:
        slab (Slab or Structure): structure or slab corresponding
            to the slab to be calculated
        transmuter (bool): whether or not to use a TransmuterFW based
            on slab params, if this option is selected, input slab must
            be a Slab object (as opposed to Structure)
        slab_gen_params (dict): dictionary of slab generation parameters
            used to generate the slab, necessary to get the slab
            that corresponds to the bulk structure
        vasp_input_set (VaspInputSet): vasp_input_set corresponding to
            the slab calculation
        parents (Fireworks or list of ints): parent FWs
        db_file (string): path to database file
        vasp_cmd (string): vasp command

    Returns:
        Firework corresponding to slab calculation
    """
    vasp_input_set = vasp_input_set or MVLSlabSet(slab)

    # If a bulk_structure is specified, generate the set of transformations,
    # else just create an optimize FW with the slab
    if transmuter:
        if not isinstance(slab, Slab):
            raise ValueError("transmuter mode requires slab to be a Slab object")

        # Get transformation from oriented bulk
        oriented_bulk = slab.oriented_unit_cell
        slab_trans_params = slab_gen_params.copy()
        slab_trans_params.update({"miller_index": [0, 0, 1],
                                  "shift": slab.shift})
        # TODO: It would be a lot more convenient if slab_gen_params could
        #       be intelligently detected from the slab itself, rather than
        #       having to be passed
        trans_struct = SlabTransformation(**slab_trans_params)
        slab_from_bulk = trans_struct.apply_transformation(oriented_bulk)

        # Ensures supercell construction
        supercell_trans = SupercellTransformation.from_scaling_factors(
            round(slab.lattice.a / slab_from_bulk.lattice.a),
            round(slab.lattice.b / slab_from_bulk.lattice.b))

        # Get adsorbates for InsertSitesTransformation
        if "adsorbate" in slab.site_properties.get("surface_properties", ""):
            ads_sites = [site for site in slab
                         if site.properties["surface_properties"] == "adsorbate"]
        else:
            ads_sites = []
        transformations = [
            "SlabTransformation", "SupercellTransformation",
            "InsertSitesTransformation", "AddSitePropertyTransformation"]
        trans_params = [slab_trans_params,
                        {"scaling_matrix": supercell_trans.scaling_matrix},
                        {"species": [site.species_string for site in ads_sites],
                         "coords": [site.frac_coords for site in ads_sites]},
                        {"site_properties": slab.site_properties}]
        fw = TransmuterFW(name=name, structure=oriented_bulk,
                          transformations=transformations,
                          transformation_params=trans_params,
                          copy_vasp_outputs=True, db_file=db_file,
                          vasp_cmd=vasp_cmd, parents=parents,
                          vasp_input_set=vasp_input_set)
    else:
        fw = OptimizeFW(name=name, structure=slab,
                        vasp_input_set=vasp_input_set, vasp_cmd=vasp_cmd,
                        db_file=db_file, parents=parents, job_type="normal")
    # Add slab metadata
    fw.tasks[-1]["additional_fields"].update({"slab": slab})
    return fw


def get_wf_slab(slab, include_bulk_opt=False, slab_gen_params=None,
                molecules=None, ads_structures_params=None, vasp_cmd="vasp",
                db_file=None, add_molecules_in_box=False):
    """
    Gets a workflow corresponding to a slab calculation along with optional
    adsorbate calcs and precursor oriented unit cell optimization

    Args:
        slabs (list of Slabs or Structures): slabs to calculate
        transmuter (bool): whether to use transmuters in constructing the
            fireworks, which will

        slab_gen_params (dict): dictionary of slab generation parameters
            used to generate the slab, necessary to get the slab
            that corresponds to the bulk structure if in that mode
        molecules (list of Molecules): molecules to place as adsorbates
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

    if molecules is None:
        molecules = []

    if ads_structures_params is None:
        ads_structures_params = {}

    if include_bulk_opt:
        oriented_bulk = slab.oriented_unit_cell
        vis = MVLSlabSet(oriented_bulk, bulk=True)
        fws.append(OptimizeFW(structure=oriented_bulk, vasp_input_set=vis,
                              vasp_cmd="vasp", db_file=db_file))
        parents = fws[-1]
    name = slab.composition.reduced_formula
    if getattr(slab, "miller_index", None):
        name += "_{}".format(slab.miller_index)
    fws.append(get_slab_fw(slab, include_bulk_opt, slab_gen_params,
                           db_file=db_file, vasp_cmd=vasp_cmd,
                           parents=parents, name=name + " slab optimization"))
    for molecule in molecules:
        ads_slabs = AdsorbateSiteFinder(slab).generate_adsorption_structures(
            molecule, **ads_structures_params)
        for n, ads_slab in enumerate(ads_slabs):
            ads_name = "{}-{} adsorbate optimization {}".format(
                molecule.composition.formula, name, n)
            fws.append(get_slab_fw(ads_slab, include_bulk_opt, slab_gen_params,
                                   db_file=db_file, vasp_cmd=vasp_cmd,
                                   parents=parents, name=ads_name))

    # TODO: add analysis framework
    if isinstance(slab, Slab):
        name = "{}_{} slab workflow".format(
            slab.composition.reduced_composition, slab.miller_index)
    else:
        name = "{} slab workflow".format(slab.composition.reduced_composition)
    wf = Workflow(fws, name=name)

    # Add optional molecules workflow
    if add_molecules_in_box:
        molecule_wf = get_wf_molecules(molecules, db_file=db_file,
                                       vasp_cmd=vasp_cmd)
        wf.append_wf(molecule_wf)
    return wf


def get_wf_molecules(molecules, vasp_input_set=None, db_file=None,
                     vasp_cmd="vasp", name=""):
    """
    Args:
        molecules (Molecules): list of molecules to calculate
        vasp_input_set (DictSet): VaspInputSet for molecules
        db_file (string): database file path
        vasp_cmd (string): VASP command
        name (string): name for workflow

    Returns:
        workflow consisting of molecule calculations
    """
    fws = []
    # TODO: this should probably include vibrations so that free energy
    #       corrections can be made
    for molecule in molecules:
        # molecule in box
        m_struct = molecule.get_boxed_structure(10, 10, 10,
                                                offset=np.array([5, 5, 5]))
        vis = vasp_input_set or MVLSlabSet(m_struct)
        fws.append(OptimizeFW(structure=molecule, job_type="normal",
                              vasp_input_set=vis, db_file=db_file,
                              vasp_cmd=vasp_cmd))
    name = name or "molecules workflow"
    return Workflow(fws, name=name)


# TODO: this will duplicate a precursor optimization for slabs with
#       the same miller index, but different shift
def get_wfs_all_slabs(bulk_structure, include_bulk_opt=False,
                      molecules=None, max_index=1, slab_gen_params=None,
                      ads_structures_params=None, vasp_cmd="vasp",
                      db_file=None, add_molecules_in_box=False):
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
        list of slab-specific Workflows
    """
    # TODO: these could be more well-thought out defaults
    sgp = slab_gen_params or {"min_slab_size": 7.0, "min_vacuum_size": 20.0}
    slabs = generate_all_slabs(bulk_structure, max_index=max_index, **sgp)
    wfs = []
    for slab in slabs:
        slab_wf = get_wf_slab(slab, include_bulk_opt, sgp, molecules,
                              ads_structures_params, vasp_cmd, db_file)
        wfs.append(slab_wf)

    if add_molecules_in_box:
        wfs.append(get_wf_molecules(molecules, db_file=db_file,
                                    vasp_cmd=vasp_cmd))

    return wfs
