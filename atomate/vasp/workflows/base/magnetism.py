# coding: utf-8

import os

from atomate.vasp.fireworks.core import OptimizeFW, StaticFW
from fireworks import Workflow, Firework
from atomate.vasp.powerups import add_tags, add_additional_fields_to_taskdocs,\
    add_wf_metadata, add_common_powerups
from atomate.vasp.workflows.base.core import get_wf
from atomate.vasp.firetasks.parse_outputs import MagneticDeformationToDB, MagneticOrderingsToDB

from pymatgen.transformations.advanced_transformations import MagOrderParameterConstraint, \
    MagOrderingTransformation
from pymatgen.analysis.local_env import MinimumDistanceNN

from atomate.utils.utils import get_logger
logger = get_logger(__name__)

from atomate.vasp.config import VASP_CMD, DB_FILE, ADD_WF_METADATA
from uuid import uuid4
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.core import Lattice, Structure
from pymatgen.analysis.magnetism.analyzer import CollinearMagneticStructureAnalyzer, Ordering

__author__ = "Matthew Horton"
__maintainer__ = "Matthew Horton"
__email__ = "mkhorton@lbl.gov"
__status__ = "Development"
__date__ = "March 2017"

__magnetic_deformation_wf_version__ = 1.2
__magnetic_ordering_wf_version__ = 1.2

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


def get_wf_magnetic_deformation(structure,
                                common_params=None,
                                vis=None):
    """
    Minimal workflow to obtain magnetic deformation proxy, as
    defined by Bocarsly et al. 2017, doi: 10.1021/acs.chemmater.6b04729

    :param structure: input structure, must be structure with magnetic
    elements, such that pymatgen will initalize ferromagnetic input by
    default -- see MPRelaxSet.yaml for list of default elements
    :param common_params (dict): Workflow config dict, in the same format
    as in presets/core.py
    :param vis: (VaspInputSet) A VaspInputSet to use for the first FW
    :return:
    """

    if not structure.is_ordered:
        raise ValueError("Please obtain an ordered approximation of the input structure.")

    structure = structure.get_primitive_structure(use_site_props=True)

    # using a uuid for book-keeping,
    # in a similar way to other workflows
    uuid = str(uuid4())

    c = {'vasp_cmd': VASP_CMD, 'db_file': DB_FILE}
    if common_params:
        c.update(common_params)

    wf = get_wf(structure, "magnetic_deformation.yaml",
                common_params=c, vis=vis)

    fw_analysis = Firework(
        MagneticDeformationToDB(
            db_file=DB_FILE,
            wf_uuid=uuid,
            to_db=c.get("to_db", True)),
        name="MagneticDeformationToDB")

    wf.append_wf(Workflow.from_Firework(fw_analysis), wf.leaf_fw_ids)

    wf = add_common_powerups(wf, c)

    if c.get("ADD_WF_METADATA", ADD_WF_METADATA):
        wf = add_wf_metadata(wf, structure)

    wf = add_additional_fields_to_taskdocs(wf, {
        'wf_meta': {
            'wf_uuid': uuid,
            'wf_name': 'magnetic_deformation',
            'wf_version': __magnetic_deformation_wf_version__
        }})

    return wf


def get_wf_magnetic_orderings(structure,
                              vasp_cmd=VASP_CMD,
                              db_file=DB_FILE,
                              default_magmoms=None,
                              attempt_ferrimagnetic=False,
                              num_orderings=10,
                              max_cell_size=None,
                              vasp_input_set_kwargs=None,
                              **kwargs):
    """
    Will try several different collinear magnetic orderings
    for a given input structure, and output a summary to a
    database or JSON file.

    If the input structure has magnetic moments defined, these
    will be used to give a hint as to which elements are
    magnetic, otherwise magnetic elements will be guessed
    (this can be changed using default_magmoms kwarg).

    Ferrimagnetic options are in beta.

    :param structure: input structure
    :param vasp_cmd: as elsewhere in atomate
    :param db_file: as elsewhere in atomate
    :param default_magmoms (dict): dict of magnetic elements
    to their initial magnetic moment in µB
    :param attempt_ferrimagnetic (bool): whether ot not to
    attempt ferrimagnetic structures
    :param num_orderings (int): This is the number of each
    type of ordering to attempt (behind the scenes, it is
    passed to pymatgen's transformation's return_ranked_list).
    Note this is per strategy, so it will return num_orderings
    AFM orderings, but also num_orderings ferrimagnetic by motif,
    etc. if attempt_ferrimagnetic == True
    :param max_cell_size (int): The max_cell_size to consider. If
    too large, enumeration will take a long time! By default will
    try a sensible value (between 4 and 1 depending on number of
    magnetic sites in primitive cell).
    :param vasp_input_set_kwargs: kwargs to pass to the
    vasp input set, default is `{'user_incar_settings':
    {'ISYM': 0, 'LASPH': True}`
    :return:
    """

    if not structure.is_ordered:
        raise ValueError("Please obtain an ordered approximation of the input structure.")

    structure = structure.get_primitive_structure(use_site_props=True)

    formula = structure.composition.reduced_formula
    msa = CollinearMagneticStructureAnalyzer(structure,
                                             default_magmoms=default_magmoms,
                                             overwrite_magmom_mode="none")

    if not msa.is_collinear:
        raise ValueError("Input structure ({}) is non-collinear.".format(formula))

    # strip out existing magmoms, can cause conflicts with transformation otherwise
    if 'magmom' in structure.site_properties:
        structure.remove_site_property('magmom')

    mag_species_spin = msa.magnetic_species_and_magmoms
    types_mag_species = msa.types_of_magnetic_specie
    types_mag_elements = {sp.symbol for sp in types_mag_species}
    num_mag_sites = msa.number_of_magnetic_sites
    num_unique_sites = msa.number_of_unique_magnetic_sites()
    max_cell_size = max_cell_size if max_cell_size else max(1, int(4/num_mag_sites))

    def _add_structures(ordered_structures, structures_to_add):
        """
        Transformations with return_ranked_list can return either
        just Structures or dicts (or sometimes lists!) -- until this
        is fixed, we use this function to concat structures given
        by the transformation.
        """
        if structures_to_add:
            # type conversion
            if isinstance(structures_to_add, Structure):
                structures_to_add = [structures_to_add]
            structures_to_add = [s["structure"] if isinstance(s, dict)
                                 else s for s in structures_to_add]
            # concatenation
            ordered_structures += structures_to_add
            logger.info('Adding {} ordered structures'.format(len(structures_to_add)))
        return ordered_structures

    # we can't enumerate all possible orderings, this combination
    # of if statements applies heuristics which work for a large
    # number of materials

    # we start with a ferromagnetic ordering
    fm_structure = msa.get_ferromagnetic_structure()
    # store magmom as spin property, to be consistent with pymatgen's transformations
    fm_structure.add_spin_by_site(fm_structure.site_properties['magmom'])
    fm_structure.remove_site_property('magmom')
    ordered_structures = [msa.get_ferromagnetic_structure()]

    # we enumerate simple AFM cases first
    constraint = MagOrderParameterConstraint(
        0.5,
        species_constraints=list(map(str, types_mag_species))
    )

    trans = MagOrderingTransformation(mag_species_spin,
                                      order_parameter=[constraint],
                                      max_cell_size=max_cell_size)
    structures_to_add = trans.apply_transformation(structure,
                                                   return_ranked_list=num_orderings)
    ordered_structures = _add_structures(ordered_structures,
                                         structures_to_add)

    # we also try ferrimagnetic orderings by motif for single magnetic species
    if attempt_ferrimagnetic and num_unique_sites > 1 and len(types_mag_elements) == 1:
        # here we try AFM, but we also try ferrimagnetic
        # orderings

        # detect co-ordination numbers
        nn = MinimumDistanceNN()
        cns = [nn.get_cn(structure, n) for n in range(len(structure))]
        is_magnetic_sites = [True if site.specie in types_mag_species
                            else False for site in structure]
        cns = [cn if is_magnetic_site else 0
               for cn, is_magnetic_site in zip(cns, is_magnetic_sites)]
        structure.add_site_property('cn', cns)

        unique_cns = set(cns) - {0}

        for cn in unique_cns:

            constraints = [
                MagOrderParameterConstraint(
                    0.5,
                    site_constraint_name='cn',
                    site_constraints=cn
                ),
                MagOrderParameterConstraint(
                    1.0,
                    site_constraint_name='cn',
                    site_constraints=list(unique_cns - {cn})
                )
            ]

            trans = MagOrderingTransformation(mag_species_spin,
                                              order_parameter=constraints,
                                              max_cell_size=max_cell_size)

            structures_to_add = trans.apply_transformation(structure,
                                                           return_ranked_list=num_orderings)

            ordered_structures = _add_structures(ordered_structures,
                                                 structures_to_add)

    # and also try ferrimagnetic when there are multiple magnetic species
    elif attempt_ferrimagnetic and len(types_mag_species) > 1:

        for sp in types_mag_species:

            constraints = [
                MagOrderParameterConstraint(
                    0.5,
                    species_constraints=str(sp)
                ),
                MagOrderParameterConstraint(
                    1.0,
                    species_constraints=list(map(str, set(types_mag_species) - {sp}))
                )
            ]

            trans = MagOrderingTransformation(mag_species_spin,
                                              order_parameter=constraints,
                                              max_cell_size=max_cell_size)

            structures_to_add = trans.apply_transformation(structure,
                                                           return_ranked_list=num_orderings)

            ordered_structures = _add_structures(ordered_structures,
                                                 structures_to_add)

    # in case we've introduced duplicates, let's remove them
    structures_to_remove = []
    for idx, ordered_structure in enumerate(ordered_structures):
        if idx not in structures_to_remove:
            matches = [ordered_structure.matches(s)
                       for s in ordered_structures]
            structures_to_remove += [match_idx for match_idx, match in enumerate(matches)
                                     if (match and idx != match_idx)]

    logger.info('Removing {} duplicate ordered structures'.format(len(structures_to_remove)))
    ordered_structures = [s for idx, s in enumerate(ordered_structures)
                          if idx not in structures_to_remove]

    # indexes keeps track of which ordering we tried first
    # it helps give us statistics for how many orderings we
    # have to try before we get the ground state
    indexes = list(range(len(ordered_structures)))

    # if our input structure isn't in our generated structures,
    # let's add it manually
    matches = [msa.matches_ordering(s) for s in ordered_structures]
    tags = []
    if not any(matches):
        ordered_structures.append(structure)
        indexes += [-1]
        logger.info("Input structure not present in enumerated structures, adding...")
    else:
        # keep a note of which structure is our input
        # this is mostly for book-keeping
        tags.append('Input structure index: {}'.format(matches.index(True)))
        logger.info("Input structure was found in enumerated "
                    "structures at index {}".format(matches.index(True)))

    fws = []
    analysis_parents = []

    for idx, ordered_structure in enumerate(ordered_structures):

        msa = CollinearMagneticStructureAnalyzer(ordered_structure)

        name = "ordering {} {} -".format(indexes[idx], msa.ordering.value)

        # get keyword arguments for VaspInputSet
        relax_vis_kwargs = {'user_incar_settings': {'ISYM': 0, 'LASPH': True}}
        if vasp_input_set_kwargs:
            relax_vis_kwargs.update(vasp_input_set_kwargs)

        if msa.ordering == Ordering.NM:
            # use with care, in general we *do not* want a non-spin-polarized calculation
            # just because initial structure has zero magnetic moments; used here for
            # calculation of magnetic deformation proxy
            relax_vis_kwargs['user_incar_settings'].update({'ISPIN': 1})

        vis = MPRelaxSet(ordered_structure, **relax_vis_kwargs)

        # relax
        fws.append(OptimizeFW(ordered_structure, vasp_input_set=vis,
                              vasp_cmd=vasp_cmd, db_file=db_file,
                              max_force_threshold=0.05,
                              half_kpts_first_relax=True,
                              name=name+" optimize"))

        # static
        fws.append(StaticFW(ordered_structure, vasp_cmd=vasp_cmd,
                            db_file=db_file,
                            name=name+" static",
                            prev_calc_loc=True, parents=fws[-1]))

        analysis_parents.append(fws[-1])

    uuid = str(uuid4())
    fw_analysis = Firework(MagneticOrderingsToDB(db_file=db_file,
                                                 wf_uuid=uuid,
                                                 auto_generated=False,
                                                 name="MagneticOrderingsToDB",
                                                 parent_structure=structure,
                                                 strategy=vasp_input_set_kwargs),
                           name="Magnetic Orderings Analysis", parents=analysis_parents)
    fws.append(fw_analysis)

    wf_name = "{} - magnetic orderings".format(formula)
    wf = Workflow(fws, name=wf_name)

    wf = add_additional_fields_to_taskdocs(wf, {
        'wf_meta': {
            'wf_uuid': uuid,
            'wf_name': 'magnetic_orderings',
            'wf_version': __magnetic_ordering_wf_version__
        }})

    tag = "magnetic_orderings group: >>{}<<".format(uuid)
    wf = add_tags(wf, [tag])

    return wf

if __name__ == "__main__":
    from fireworks import LaunchPad

    latt = Lattice.cubic(4.17)
    species = ["Ni", "O"]
    coords = [[0.00000, 0.00000, 0.00000],
              [0.50000, 0.50000, 0.50000]]
    NiO = Structure.from_spacegroup(225, latt, species, coords)

    wf_deformation = get_wf_magnetic_deformation(NiO)

    wf_orderings = get_wf_magnetic_orderings(NiO)

    lpad = LaunchPad.auto_load()
    lpad.add_wf(wf_orderings)
    lpad.add_wf(wf_deformation)
