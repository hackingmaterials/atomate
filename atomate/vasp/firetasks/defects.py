# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module provides ability to setup and parse defects on fly...
Note alot of this structure is taken from pycdt.utils.parse_calculations and pycdt.......

Requirements:
 - bulk calculation is finished (and vasprun.xml + Locpot)
 - a dielectric constant/tensor is provided
 - defect+chg calculation is finished 

 Soft requirements:
 	- Bulk and defect OUTCAR files (if charge correction by Kumagai et al. is desired)
 	- Hybrid bulk bandstructure / simple bulk structure calculation (if bandshifting is desired)
"""

import os
import itertools
import numpy as np

from monty.json import jsanitize

from pymatgen.io.vasp import Vasprun, Locpot, Poscar
from pymatgen import MPRester
from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPRelaxSet, MVLScanRelaxSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.analysis.defects.generators import VacancyGenerator, SubstitutionGenerator, \
    InterstitialGenerator, VoronoiInterstitialGenerator, SimpleChargeGenerator

from fireworks import FiretaskBase, FWAction, explicit_serialize

from atomate.utils.utils import get_logger
from atomate.vasp.fireworks.core import TransmuterFW

from monty.serialization import dumpfn
from monty.json import MontyEncoder

logger = get_logger(__name__)

def optimize_structure_sc_scale(inp_struct, final_site_no):
    """
    A function for finding optimal supercell transformation
    by maximizing the nearest image distance with the number of
    atoms remaining less than final_site_no

    Args:
        inp_struct: input pymatgen Structure object
        final_site_no (float or int): maximum number of atoms
        for final supercell

    Returns:
        3 x 1 array for supercell transformation
    """
    if final_site_no <= len(inp_struct.sites):
        final_site_no = len(inp_struct.sites)

    dictio={}
    #consider up to a 7x7x7 supercell
    for kset in itertools.product(range(1,7), range(1,7), range(1,7)):
        num_sites = len(inp_struct) * np.product(kset)
        if num_sites > final_site_no:
            continue

        struct = inp_struct.copy()
        struct.make_supercell(kset)

        #find closest image
        min_dist = 1000.
        for image_array in itertools.product( range(-1,2), range(-1,2), range(-1,2)):
            if image_array == (0,0,0):
                continue
            distance = struct.get_distance(0, 0, image_array)
            if distance < min_dist:
                min_dist = distance

        min_dist = round(min_dist, 3)
        if min_dist in dictio.keys():
            if dictio[min_dist]['num_sites'] > num_sites:
                dictio[min_dist].update( {'num_sites': num_sites, 'supercell': kset[:]})
        else:
            dictio[min_dist] = {'num_sites': num_sites, 'supercell': kset[:]}

    if not len(dictio.keys()):
        raise RuntimeError('could not find any supercell scaling vector')

    min_dist = max( list(dictio.keys()))
    biggest = dictio[ min_dist]['supercell']

    return biggest


@explicit_serialize
class DefectSetupFiretask(FiretaskBase):
    """
    Run defect supercell setup

    Args:
        structure (Structure): input structure to have defects run on
        cellmax (int): maximum supercell size to consider for supercells
        conventional (bool):
            flag to use conventional structure (rather than primitive) for supercells,
            defaults to True.
        vasp_cmd (string):
            the vasp cmd
        db_file (string):
            the db file
        user_incar_settings (dict):
            a dictionary of incar settings specified by user for both bulk and defect supercells
            note that charges do not need to be set in this dicitionary
        user_kpoints_settings (dict or Kpoints pmg object):
            a dictionary of kpoint settings specific by user OR an Actual Kpoint set to be used for the calculation

        vacancies (list):
            If list is totally empty, all vacancies are considered (default).
            If only specific vacancies are desired then add desired Element symbol to the list
                ex. ['Ga'] in GaAs structure will only produce Galium vacancies

            if NO vacancies are desired, then just add an empty list to the list
                ex. [ [] ]  yields no vacancies

        substitutions (dict):
            If dict is totally empty, all intrinsic antisites are considered (default).
            If only specific antisites/substituions are desired then add vacant site type as key, with list of
                sub site symbol as value
                    ex 1. {'Ga': ['As'] } in GaAs structure will only produce Arsenic_on_Gallium antisites
                    ex 2. {'Ga': ['Sb'] } in GaAs structure will only produce Antimonide_on_Gallium substitutions

            if NO antisites or substitutions are desired, then just add an empty dict
                ex. {'None':{}}  yields no antisites or subs

        interstitials (list):
            If list is totally empty, NO interstitial defects are considered (default).
            Option 1 for generation: If one wants to use Pymatgen to predict interstitial
                    then list of pairs of [symbol, generation method (str)] can be provided
                        ex. ['Ga', 'Voronoi'] in GaAs structure will produce Galium interstitials from the
                            Voronoi site finding algorithm
                        NOTE: only options for interstitial generation are "Voronoi" and "InFit"
            Option 2 for generation: If user wants to add their own interstitial sites for consideration
                    the list of pairs of [symbol, Interstitial object] can be provided, where the
                    Interstitial pymatgen.analysis.defects.core object is used to describe the defect of interest
                    NOTE: use great caution with this approach. You better be sure that the supercell with Interstitial in it
                        is same as the bulk supercell...

        initial_charges (dict):
            says how to specify initial charges for each defect.
            An empty dict (DEFAULT) is to do a fairly restrictive charge generation method:
                for vacancies: use bond valence method to assign oxidation states and consider
                    negative of the vacant site's oxidation state as single charge to try
                antisites and subs: use bond valence method to assign oxidation states and consider
                    negative of the vacant site's oxidation state as single charge to try +
                    added to likely charge of substitutional site (closest to zero)
                interstitial: charge zero
            For non empty dict, charges are specified as:
                initial_charges = {'vacancies': {'Ga': [-3,2,1,0]},
                                   'substitutions': {'Ga': {'As': [0]} },
                                   'interstitials': {}}
                in the GaAs structure this makes vacancy charges in states -3,-2,-1,0; Ga_As antisites in the q=0 state,
                and all other defects will have charges generated in the restrictive automated format stated for DEFAULT

    """
    def run_task(self, fw_spec):
        if os.path.exists("POSCAR"):
            structure =  Poscar.from_file("POSCAR").structure
        else:
            structure = self.get("structure")

        if self.get("conventional", True):
            structure = SpacegroupAnalyzer(structure).get_conventional_standard_structure()

        fws, parents = [], []

        cellmax=self.get("cellmax", 128)
        sc_scale = optimize_structure_sc_scale(structure, cellmax)

        #First Firework is for bulk supercell
        bulk_supercell = structure.copy()
        bulk_supercell.make_supercell(sc_scale)
        num_atoms = len(bulk_supercell)

        user_incar_settings = self.get("user_incar_settings", {})
        user_kpoints_settings = self.get("user_kpoints_settings", {})

        bulk_incar_settings = {"EDIFF":.0001, "EDIFFG": 0.001, "ISMEAR":0, "SIGMA":0.05, "NSW": 0, "ISIF": 2,
                               "ISPIN":2,  "ISYM":2, "LVHAR":True, "LVTOT":True, "LWAVE": True}
        bulk_incar_settings.update( user_incar_settings)

        reciprocal_density = 100
        kpoints_settings = user_kpoints_settings if user_kpoints_settings else {"reciprocal_density": reciprocal_density}
        vis = MPRelaxSet( bulk_supercell,
                          user_incar_settings=bulk_incar_settings,
                          user_kpoints_settings=kpoints_settings)

        supercell_size = sc_scale * np.identity(3)
        bulk_tag = "{}:bulk_supercell_{}atoms".format(structure.composition.reduced_formula, num_atoms)
        stat_fw = TransmuterFW(name = bulk_tag, structure=structure,
                               transformations=['SupercellTransformation'],
                               transformation_params=[{"scaling_matrix": supercell_size}],
                               vasp_input_set=vis, copy_vasp_outputs=False,
                               vasp_cmd=self.get("vasp_cmd", ">>vasp_cmd<<"),
                               db_file=self.get("db_file", ">>db_file<<"))

        fws.append(stat_fw)

        # make defect set
        vacancies = self.get("vacancies", list())
        substitutions = self.get("substitutions", dict())
        interstitials = self.get("interstitials", list())
        initial_charges  = self.get("initial_charges", dict())

        def_structs = []
        #a list with following dict structure for each entry:
        # {'defect': pymatgen defect object type,
        # 'charges': list of charges to run}

        if not vacancies:
            #default: generate all vacancies...
            b_struct = structure.copy()
            VG = VacancyGenerator( b_struct)
            for vac_ind, vac in enumerate(VG):
                vac_symbol = vac.site.specie.symbol

                charges = []
                if initial_charges:
                    if 'vacancies' in initial_charges.keys():
                        if vac_symbol in initial_charges['vacancies']:
                            #NOTE if more than one type of vacancy for a given specie, this will assign same charges to all
                            charges = initial_charges['vacancies'][vac_symbol]

                if not len(charges):
                    SCG = SimpleChargeGenerator(vac.copy())
                    charges = [v.charge for v in SCG]

                def_structs.append({'charges': charges, 'defect': vac.copy()})

        else:
            #only create vacancies of interest...
            for elt_type in vacancies:
                b_struct = structure.copy()
                VG = VacancyGenerator( b_struct)
                for vac_ind, vac in enumerate(VG):
                    vac_symbol = vac.site.specie.symbol
                    if elt_type != vac_symbol:
                        continue

                    charges = []
                    if initial_charges:
                        if 'vacancies' in initial_charges.keys():
                            if vac_symbol in initial_charges['vacancies']:
                                #NOTE if more than one type of vacancy for a given specie, this will assign same charges to all
                                charges = initial_charges['vacancies'][vac_symbol]

                    if not len(charges):
                        SCG = SimpleChargeGenerator(vac.copy())
                        charges = [v.charge for v in SCG]

                    def_structs.append({'charges': charges, 'defect': vac.copy()})


        if not substitutions:
            #default: set up all intrinsic antisites
            for sub_symbol in [elt.symbol for elt in bulk_supercell.types_of_specie]:
                b_struct = structure.copy()
                SG = SubstitutionGenerator(b_struct, sub_symbol)
                for as_ind, sub in enumerate(SG):
                    #find vac_symbol to correctly label defect
                    poss_deflist = sorted(sub.bulk_structure.get_sites_in_sphere(sub.site.coords, 2, include_index=True), key=lambda x: x[1])
                    defindex = poss_deflist[0][2]
                    vac_symbol = sub.bulk_structure[defindex].specie.symbol

                    charges = []
                    if initial_charges:
                        if 'substitutions' in initial_charges.keys():
                            if vac_symbol in initial_charges['substitutions']:
                                #NOTE if more than one type of substituion for a given specie, this will assign same charges to all
                                if sub_symbol in initial_charges['substitutions'][vac_symbol].keys():
                                    charges = initial_charges['substitutions'][vac_symbol][sub_symbol]
                    if not len(charges):
                        SCG = SimpleChargeGenerator(sub.copy())
                        charges = [v.charge for v in SCG]

                    def_structs.append({'charges': charges, 'defect': sub.copy()})
        else:
            #only set up specified antisite / substituion types
            for vac_symbol, sub_list in substitutions.items():
                for sub_symbol in sub_list:
                    b_struct = structure.copy()
                    SG = SubstitutionGenerator(b_struct, sub_symbol)
                    for as_ind, sub in enumerate(SG):
                        #find vac_symbol for this sub defect
                        poss_deflist = sorted(sub.bulk_structure.get_sites_in_sphere(sub.site.coords, 2, include_index=True), key=lambda x: x[1])
                        defindex = poss_deflist[0][2]
                        gen_vac_symbol = sub.bulk_structure[defindex].specie.symbol
                        if vac_symbol != gen_vac_symbol: #only consider subs on specfied vac_symbol site
                            continue

                        charges = []
                        if initial_charges:
                            if 'substitutions' in initial_charges.keys():
                                if vac_symbol in initial_charges['substitutions']:
                                    #NOTE if more than one type of substituion for a given specie, this will assign same charges to all
                                    if sub_symbol in initial_charges['substitutions'][vac_symbol].keys():
                                        charges = initial_charges['substitutions'][vac_symbol][sub_symbol]
                        if not len(charges):
                            SCG = SimpleChargeGenerator(sub.copy())
                            charges = [v.charge for v in SCG]

                        def_structs.append({'charges': charges, 'defect': sub.copy()})


        if interstitials:
            #default = do not include interstitial defects

            def get_charges_from_inter( inter_elt):
                inter_charges = []
                if initial_charges:
                    if 'interstitials' in initial_charges.keys():
                        if inter_elt in initial_charges['interstitials']:
                            #NOTE if more than one type of interstitial for a given specie, this will assign same charges to all
                            inter_charges = initial_charges['interstitials'][inter_elt]

                if not len(inter_charges):
                    SCG = SimpleChargeGenerator(inter_elt)
                    inter_charges = [v.charge for v in SCG]
                return inter_charges

            for elt_type, elt_val in interstitials:
                if type(elt_val) == str:
                    b_struct = structure.copy()
                    if elt_val == 'Voronoi':
                        IG = VoronoiInterstitialGenerator(b_struct, elt_type)
                    elif elt_val == 'InFit':
                        IG = InterstitialGenerator(b_struct, elt_type)
                    else:
                        raise ValueError('Interstitial finding method not recognized. '
                                         'Please choose either Voronoi or InFit.')

                    for inter_ind, inter in enumerate(IG):
                        charges = get_charges_from_inter( inter)
                        def_structs.append({'charges': charges, 'defect': inter.copy()})
                else:
                    charges = get_charges_from_inter( elt_type)
                    def_structs.append({'charges': charges, 'defect': elt_val.copy()})

        stdrd_defect_incar_settings = {"EDIFF": 0.0001, "EDIFFG": 0.001, "IBRION":2, "ISMEAR":0, "SIGMA":0.05,
                                       "ISPIN":2,  "ISYM":2, "LVHAR":True, "LVTOT":True, "NSW": 100,
                                       "NELM": 60, "ISIF": 2, "LAECHG":False, "LWAVE": True}
        stdrd_defect_incar_settings.update( user_incar_settings)

        # now that def_structs is assembled, set up Transformation FW for all defect + charge combinations
        for defcalc in def_structs:
            #get defect supercell and defect site for parsing purposes
            defect = defcalc['defect'].copy()
            defect_sc = defect.generate_defect_structure( supercell = supercell_size)
            struct_for_defect_site = Structure(defect.bulk_structure.copy().lattice,
                                               [defect.site.specie],
                                               [defect.site.frac_coords],
                                               to_unit_cell=True, coords_are_cartesian=False)
            struct_for_defect_site.make_supercell(supercell_size)
            defect_site = struct_for_defect_site[0]

            #iterate over all charges to be run
            for charge in defcalc['charges']:
                chgdstruct = defect_sc.copy()
                chgdstruct.set_charge(charge)  #NOTE that the charge will be reflected in NELECT of INCAR because use_structure_charge=True

                reciprocal_density = 100
                kpoints_settings = user_kpoints_settings if user_kpoints_settings else {"reciprocal_density": reciprocal_density}
                defect_input_set = MPRelaxSet( chgdstruct,
                                               user_incar_settings=stdrd_defect_incar_settings.copy(),
                                               user_kpoints_settings=kpoints_settings,
                                               use_structure_charge=True)

                defect_for_trans_param = defect.copy()
                defect_for_trans_param.set_charge(charge)
                chgdef_trans = ["DefectTransformation"]
                chgdef_trans_params = [{"scaling_matrix": supercell_size,
                                        "defect": defect_for_trans_param}]

                def_tag = "{}:{}_{}_{}atoms".format(structure.composition.reduced_formula,
                                                      defect.name, charge, num_atoms)
                fw = TransmuterFW( name = def_tag, structure=structure,
                                   transformations=chgdef_trans,
                                   transformation_params=chgdef_trans_params,
                                   vasp_input_set=defect_input_set,
                                   vasp_cmd=self.get("vasp_cmd", ">>vasp_cmd<<"),
                                   copy_vasp_outputs=False,
                                   db_file=self.get("db_file", ">>db_file<<"),
                                   bandstructure_mode="auto",
                                   defect_wf_parsing=defect_site)

                fws.append(fw)

        return FWAction(detours=fws)
