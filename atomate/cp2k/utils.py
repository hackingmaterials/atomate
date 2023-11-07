from pymatgen.analysis.defects.generators import VacancyGenerator, InterstitialGenerator, SubstitutionGenerator
from pymatgen.analysis.defects.core import Vacancy
from pymatgen.core.structure import Structure
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
        for _vac in VG:
            vac = GhostVacancy(_vac.bulk_structure, _vac.site, _vac.charge, _vac.multiplicity)
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

    # TODO shouldn't this be done in defect object?
    for d in defect_structures:
        d.bulk_structure.set_charge(d.bulk_structure.charge + d.charge)

    return defect_structures


def optimize_structure_sc_scale(inp_struct, final_site_no=None):
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


def optimize_structure_sc_scale_by_length(inp_struct, minimum_distance):
    """
    Given a structure, provide the scaling matrix such that the points in the lattice maintain
    a minimum distance from one another.

    Args:
        inp_struct (Structure): structure to consider
        minimum_distance (float or int): minimum distance desired between atoms. Note that if the
            value of minimum distance is < 0, this means a scalar will be returned, defining the
            minimum uniform scaling required to achieve minimum_distance. If minimum distance is
            > 0, then a vector will be returned, giving the appropriate scaling in each direction
            to achieve minimum_distance.
    :return:
    """
    if minimum_distance == 0:
        raise ValueError("Cannot assign a minimum distance of 0!")
    elif minimum_distance < 0:
        minimum_distance = -minimum_distance
        return np.ceil(minimum_distance/np.linalg.norm(inp_struct.lattice.matrix, axis=0)).astype(int)
    else:
        return np.ceil(minimum_distance/np.linalg.norm(inp_struct.lattice.matrix, axis=0)).astype(int)


class GhostVacancy(Vacancy):
    """
    Current workaround for the Vacancy class in CP2K. Vacancies are normally just structures
    with an atom removed, but with CP2K we want to retain the site and turn off its interaction
    potential (Ghost atom) in order to avoid Basis set superposition error for localized basis.
    """

    def generate_defect_structure(self, supercell=(1, 1, 1)):
        """
        Returns Defective Vacancy structure, decorated with charge
        Args:
            supercell (int, [3x1], or [[]] (3x3)): supercell integer, vector, or scaling matrix
        """
        defect_structure = self.bulk_structure.copy()
        defect_structure.make_supercell(supercell)

        # create a trivial defect structure to find where supercell transformation moves the lattice
        struct_for_defect_site = Structure(self.bulk_structure.copy().lattice,
                                           [self.site.specie],
                                           [self.site.frac_coords],
                                           to_unit_cell=True)
        struct_for_defect_site.make_supercell(supercell)
        defect_site = struct_for_defect_site[0]

        poss_deflist = sorted(
            defect_structure.get_sites_in_sphere(defect_site.coords, 0.1, include_index=True), key=lambda x: x[1])
        defindex = poss_deflist[0][2]
        defect_structure.add_site_property('ghost', [True if i == defindex else False for i in range(len(defect_structure))])
        defect_structure.set_charge(self.charge)
        return defect_structure


