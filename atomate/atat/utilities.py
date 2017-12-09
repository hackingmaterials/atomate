# coding: utf-8

from __future__ import absolute_import, division, print_function, \
    unicode_literals

from pymatgen.core import Structure
from pymatgen.transformations.standard_transformations import DiscretizeOccupanciesTransformation
from pymatgen.analysis.structure_matcher import StructureMatcher, OccupancyComparator

from atomate.atat.database import McsqsCalcDb

__author__ = 'Matthew Horton'
__credits__ = 'Josua Vieten'
__email__ = 'mkhorton@lbl.gov'


def preprocess_disordered_structure(disordered_struct, max_denominator=8):
    """
    Sanitizes input structure to a common format: makes primitive,
    normalizes volume to unity, and rescales occupancies.

    :param disordered_struct: disordered Structure
    :param max_denominator: as in pymatgen's
    DiscretizeOccupancyTransformation
    :return: disordered Structure
    """

    # TODO: substitute elements for anonymous counterparts ?

    if disordered_struct.is_ordered:
        raise ValueError("Please only use with disordered structures.")

    disordered_struct = disordered_struct.get_primitive_structure()

    disordered_struct.scale_lattice(1.0)

    if max_denominator:
        trans = DiscretizeOccupanciesTransformation(max_denominator)
        disordered_struct = trans.apply_transformation(disordered_struct)

    return disordered_struct


def retrieve_bestsqs_template(disordered_struct, db_file,
                              size=None, clusters=None, user_input_settings=None,
                              max_denominator=8):
    """
    Will check SQS database for a matching structure, and if found
    will return a list of bestsqs template structures along
    with associated metadata.

    Post-processing is required to map elements in the
    input `disordered_struct` to the elements in the retrieved
    `bestsqs` template.

    Returned list is of dictionaries with keys: 'disordered',
    'bestsqs', 'size', 'clusters', 'user_input_settings'.

    :param disordered_struct: disordered Structure
    :param db_file: for db containing 'sqs' collection
    :param size: if specific size required
    :param clusters: if specific clusters required
    :param user_input_settings: if specific user inputs required
    :param max_denominator: used to sanitize input structure,
    discretizes input occupancies
    :return:
    """

    disordered_struct = preprocess_disordered_structure(disordered_struct,
                                                        max_denominator=max_denominator)

    # query database
    query = {
        'spacegroup': disordered_struct.get_space_group_info()[0],
        'input.n_disordered': len(disordered_struct)
    }

    # if user specifies non-default options, make sure we check those too
    if size:
        query['size'] = size
    if clusters:
        query['input.clusters'] = clusters
    if user_input_settings:
        query['input.user_input_settings'] = user_input_settings

    db = McsqsCalcDb.from_db_file(db_file=db_file)

    docs = list(db.collection.find(query, [
        'size',
        'input.disordered',
        'input.user_input_settings',
        'input.clusters',
        'output.bestsqs',
        'output.objective_function'
    ]))

    sm = StructureMatcher(scale=True, comparator=OccupancyComparator())
    structures = [Structure.from_dict(doc['input']['disordered']) for doc in docs]
    matches = [sm.fit(disordered_struct, s) for s in structures]

    if any(matches):

        bestsqs_templates = []

        for match_index, match in enumerate(matches):
            if match:

                disordered = structures[match_index],
                bestsqs = docs[match_index]['output']['bestsqs']
                size = docs[match_index]['size']
                clusters = docs[match_index]['input']['clusters']
                user_input_settings = docs[match_index]['input']['user_input_settings']
                objective_function = docs[match_index]['output']['objective_function']

                bestsqs_templates.append({
                    'disordered': disordered,
                    'bestsqs': Structure.from_dict(bestsqs),
                    'size': size,
                    'clusters': clusters,
                    'user_input_settings': user_input_settings,
                    'objective_function': objective_function
                })

        return bestsqs_templates

    else:

        return None