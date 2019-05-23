def add_fix_two_atom_selective_dynamics(structure, fixed_index=0):
    """
        Returns structure with selective dynamics assigned to fix the position of two sites.
        Two sites will be fixed: 1) the site specified by fixed_index and 2) the site positioned furthest from the specified fixed_index site.

        Args:
            structure (Structure): Input structure (e.g. host lattice with one working ion intercalated)
            fixed_index (int): Index of site in structure whose position will be fixed (e.g. working ion site)
        Returns:
            Structure
    """
    sd_structure = structure.copy()
    sd_array = [[True, True, True]] * sd_structure.num_sites
    sd_array[fixed_index] = [False, False, False]
    ref_site = sd_structure.sites[fixed_index]
    distances = []
    for site in sd_structure.sites:
        distances.append(site.distance(ref_site))
    farthest_index = distances.index(max(distances))
    sd_array[farthest_index] = [False, False, False]
    sd_structure.add_site_property('selective_dynamics', sd_array)
    return sd_structure