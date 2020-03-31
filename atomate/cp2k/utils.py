from pymatgen.analysis.defects.generators import VacancyGenerator, InterstitialGenerator, SubstitutionGenerator
import itertools
import numpy as np

def get_defect_structures(structure, defect_dict):

    vacancies = defect_dict.get('vacancies')
    substitutions = defect_dict.get('substitutions')
    interstitials = defect_dict.get('interstititals')
    defect_structures = []

    if vacancies:
        # only create vacancies of interest...
        b_struct = structure.copy()  # base structure (un-defective)
        VG = VacancyGenerator(b_struct, include_bv_charge=False)
        for vac in VG:
            vac_ind = vac.site.specie.symbol
            if vac_ind not in vacancies.keys():
                continue
            else:
                if isinstance(vacancies[vac_ind], list):
                    for c in vacancies[vac_ind]:
                        v = vac.copy()
                        v.set_charge(c)
                        defect_structures.append(v)
                else:
                    defect_structures.append(vac)

    return defect_structures


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


from pymatgen.core.structure import Structure
from pymatgen.ext.matproj import MPRester

with MPRester() as mp:
    struc = mp.get_structure_by_material_id('mp-149')
