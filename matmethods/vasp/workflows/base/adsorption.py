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

from pymatgen.analysis.adsorption import generate_decorated_slabs,\
        AdsorbateSiteFinder
from pymatgen.transformations.advanced_transformations import SlabTransformation
from pymatgen.transformations.standard_transformations import SupercellTransformation
from pymatgen.transformations.site_transformations import InsertSitesTransformation
from pymatgen.io.vasp.sets import MVLSlabSet
from pymatgen import Structure

__author__ = 'Joseph Montoya'
__email__ = 'montoyjh@lbl.gov'

logger = get_logger(__name__)

@explicit_serialize
class PassSlabEnergy(FireTaskBase):
    """
    Placeholder just in case I need to pass something
    """

    def run_task(self, fw_spec):
        pass
        return FWAction()


@explicit_serialize
class AnalyzeAdsorption(FireTaskBase):
    """
    Analyzes the adsorption energies in a workflow
    """

    required_params = ['structure']
    optional_params = ['db_file']

    def run_task(self, fw_spec):
        pass
        """
        # Get optimized structure
        # TODO: will this find the correct path if the workflow is rerun from the start?
        optimize_loc = fw_spec["calc_locs"][0]["path"]
        logger.info("PARSING INITIAL OPTIMIZATION DIRECTORY: {}".format(optimize_loc))
        drone = VaspDrone()
        optimize_doc = drone.assimilate(optimize_loc)
        opt_struct = Structure.from_dict(optimize_doc["calcs_reversed"][0]["output"]["structure"])
        
        d = {"analysis": {}, "deformation_tasks": fw_spec["deformation_tasks"],
             "initial_structure": self['structure'].as_dict(), 
             "optimized_structure": opt_struct.as_dict()}

        # Save analysis results in json or db
        db_file = env_chk(self.get('db_file'), fw_spec)
        if not db_file:
            with open("adsorption.json", "w") as f:
                f.write(json.dumps(d, default=DATETIME_HANDLER))
        else:
            db = MMDb.from_db_file(db_file, admin=True)
            db.collection = db.db["adsorption"]
            db.collection.insert_one(d)
            logger.info("ADSORPTION ANALYSIS COMPLETE")
        return FWAction()
        """


def get_wf_adsorption(structure, adsorbate_config, vasp_input_set=None,
                      min_slab_size=7.0, min_vacuum_size=12.0,
                      max_normal_search=1, center_slab=True,
                      slab_gen_config=None, vasp_cmd="vasp",
                      slab_input_params=None, adsorbate_input_params=None,
                      db_file=None):
    """
    Returns a workflow to calculate adsorption structures and surfaces.

    Firework 1 : write vasp input set for structural relaxation,
                 run vasp,
                 pass run location,
                 database insertion.

    Firework 2 - N: Optimize Slab and Adsorbate Structures
    
    Args:
        structure (Structure): input structure to be optimized and run
        # TODO: rethink configuration
        adsorption_config_dict (dict): configuration dictionary for adsorption
            should be a set of keys corresponding to miller indices/molecules
            e.g. {"111":Molecule("CO", [[0, 0, 0], [0, 0, 1]])}
        vasp_input_set (DictVaspInputSet): vasp input set.
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.

    Returns:
        Workflow
    """

    v = vasp_input_set or MVLSlabSet(structure, bulk=True)
    fws = []
    fws.append(OptimizeFW(structure=structure, vasp_input_set=v,
                          vasp_cmd=vasp_cmd, db_file=db_file))

    max_index = max([int(i) for i in ''.join(adsorbate_config.keys())])
    slabs = generate_decorated_slabs(structure, max_index, min_slab_size,
                                     min_vacuum_size, max_normal_search,
                                     center_slab)
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
            trans = [SlabTransformation(slab.miller_index, min_slab_size, 
                                        min_vacuum_size, slab.shift, 
                                        slab_gen_config)]
            # TODO: name these more intelligently
            fws.append(TransmuterFW(name="slab calculation",
                                    structure = structure,
                                    transformations = trans,
                                    copy_vasp_outputs=True,
                                    db_file=db_file,
                                    vasp_cmd=vasp_cmd,
                                    parents=fws[0],
                                    vasp_input_set = "MVLSlabSet",
                                    vasp_input_params = slab_input_params)
                      )
            # Generate adsorbate configurations and add fws to workflow
            asf = AdsorbateSiteFinder(slab, selective_dynamics=True)
            for molecule in adsorbate_config[mi_string]:
                structures = asf.generate_adsorption_structures(molecule)

                for struct in structures:
                    trans_slab = [SlabTransformation(struct.miller_index, min_slab_size, 
                                        min_vacuum_size, slab.shift, 
                                        slab_gen_config)]
                    trans_supercell = [SupercellTransformation.from_scaling_factors(
                        struct.lattice.a / slab.lattice.a, struct.lattice.b / slab.lattice.b)]
                    ads_sites = [site for site in slab if 
                                 site.properties["surface_properties"]=="adsorbate"]
                    trans_add_adsorbate = InsertSitesTransformation(
                        [site.species_string for site in ads_sites],
                        [site.frac_coords for site in ads_sites])

                    ads_trans = [trans_slab, trans_supercell, trans_add_adsorbate]
                    # Might need to generate adsorbate input set
                    fws.append(TransmuterFW(name="adsorbate calculation",
                                    structure = structure,
                                    transformations = ads_trans,
                                    copy_vasp_outputs=True,
                                    db_file=db_file,
                                    vasp_cmd=vasp_cmd,
                                    parents=fws[0],
                                    vasp_input_set = "MVLSlabSet",
                                    vasp_input_params = adsorbate_input_params))

    wfname = "{}:{}".format(structure.composition.reduced_formula, "Adsorbate calculations")
    return Workflow(fws, name=wfname)

if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest
    from pymatgen import Molecule
    co = Molecule("CO", [[0, 0, 0], [0, 0, 1.9]])
    adsorbate_config = {"111":co}
    structure = PymatgenTest.get_structure("Si")
    wf = get_wf_adsorption(structure, adsorbate_config)
