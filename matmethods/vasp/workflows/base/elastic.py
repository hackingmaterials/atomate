# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the elastic workflow
"""

import json
from decimal import Decimal

import numpy as np

from fireworks import FireTaskBase, Firework, FWAction, Workflow
from fireworks.utilities.fw_serializers import DATETIME_HANDLER
from fireworks.utilities.fw_utilities import explicit_serialize

from matmethods.utils.utils import env_chk, get_logger
from matmethods.vasp.drones import VaspDrone
from matmethods.vasp.database import MMDb
from matmethods.vasp.fireworks.core import OptimizeFW, TransmuterFW

from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.analysis.elasticity.strain import Deformation, IndependentStrain
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.analysis.elasticity import reverse_voigt_map
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.sets import MPRelaxSet, DictSet
from pymatgen import Structure

__author__ = 'Shyam Dwaraknath, Joseph Montoya'
__email__ = 'shyamd@lbl.gov, montoyjh@lbl.gov'

logger = get_logger(__name__)

@explicit_serialize
class PassStressStrainData(FireTaskBase):
    """
    Passes the stress and deformation for an elastic deformation calculation

    Required params:
        deformation: the deformation gradient used in the elastic analysis.
    """

    required_params = ["deformation"]

    def run_task(self, fw_spec):
        v = Vasprun('vasprun.xml.gz')
        stress = v.ionic_steps[-1]['stress']
        defo = self['deformation']
        d_ind = np.nonzero(defo - np.eye(3))
        delta = Decimal((defo - np.eye(3))[d_ind][0])
        # Shorthand is d_X_V, X is voigt index, V is value
        dtype = "_".join(["d", str(reverse_voigt_map[d_ind][0]),
                          "{:.0e}".format(delta)])
        strain = IndependentStrain(defo)
        defo_dict = {'deformation_matrix': defo,
                     'strain': strain.tolist(),
                     'stress': stress}

        return FWAction(mod_spec=[{'_set': {
            'deformation_tasks->{}'.format(dtype): defo_dict}}])


@explicit_serialize
class AnalyzeStressStrainData(FireTaskBase):
    """
    Analyzes the stress/strain data of an elastic workflow to produce
    an elastic tensor and various other quantities.
    """

    required_params = ['structure']
    optional_params = ['db_file']

    def run_task(self, fw_spec):

        # Get optimized structure
        # TODO: will this find the correct path if the workflow is rerun from the start?
        optimize_loc = fw_spec["calc_locs"][0]["path"]
        logger.info("PARSING INITIAL OPTIMIZATION DIRECTORY: {}".format(optimize_loc))
        drone = VaspDrone()
        optimize_doc = drone.assimilate(optimize_loc)
        opt_struct = Structure.from_dict(
            optimize_doc["calcs_reversed"][0]["output"]["structure"])
        
        d = {"analysis": {}, "deformation_tasks": fw_spec["deformation_tasks"],
             "initial_structure": self['structure'].as_dict(), 
             "optimized_structure": opt_struct.as_dict()}

        dtypes = fw_spec["deformation_tasks"].keys()
        defos = [fw_spec["deformation_tasks"][dtype]["deformation_matrix"]
                 for dtype in dtypes]
        stresses = [fw_spec["deformation_tasks"][dtype]["stress"] for dtype in dtypes]
        stress_dict = {IndependentStrain(defo) : Stress(stress) for defo, stress 
                       in zip(defos, stresses)}
        
        logger.info("ANALYZING STRESS/STRAIN DATA")
        # DETERMINE IF WE HAVE 6 "UNIQUE" deformations
        if len(set([de[:3] for de in dtypes])) == 6:
            # Perform Elastic tensor fitting and analysis
            result = ElasticTensor.from_stress_dict(stress_dict)
            d["elastic_tensor"] = result.voigt.tolist()
            kg_average = result.kg_average
            d.update({"K_Voigt": kg_average[0], "G_Voigt": kg_average[1],
                      "K_Reuss": kg_average[2], "G_Reuss": kg_average[3],
                      "K_Voigt_Reuss_Hill": kg_average[4],
                      "G_Voigt_Reuss_Hill": kg_average[5]})
            d["universal_anisotropy"] = result.universal_anisotropy
            d["homogeneous_poisson"] = result.homogeneous_poisson

        else:
            raise ValueError("Fewer than 6 unique deformations")

        d["state"] = "successful"

        # Save analysis results in json or db
        db_file = env_chk(self.get('db_file'), fw_spec)
        if not db_file:
            with open("elasticity.json", "w") as f:
                f.write(json.dumps(d, default=DATETIME_HANDLER))
        else:
            db = MMDb.from_db_file(db_file, admin=True)
            db.collection = db.db["elasticity"]
            db.collection.insert_one(d)
            logger.info("ELASTIC ANALYSIS COMPLETE")

        return FWAction()

def get_wf_elastic_constant(structure, vasp_input_set=None, vasp_cmd="vasp", 
                            norm_deformations=[-0.01, -0.005, 0.005, 0.01],
                            shear_deformations=[-0.06, -0.03, 0.03, 0.06],
                            db_file=None, reciprocal_density=None):
    """
    Returns a workflow to calculate elastic constants.

    Firework 1 : write vasp input set for structural relaxation,
                 run vasp,
                 pass run location,
                 database insertion.

    Firework 2 - 25: Optimize Deformed Structure
    
    Firework 26: Analyze Stress/Strain data and fit the elastic tensor

    Args:
        structure (Structure): input structure to be optimized and run
        norm_deformations (list): list of values to for normal deformations
        shear_deformations (list): list of values to for shear deformations
        vasp_input_set (DictVaspInputSet): vasp input set.
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.
        reciprocal_density (int): k-points per reciprocal atom by volume

    Returns:
        Workflow
    """

    v = vasp_input_set or MPRelaxSet(structure, force_gamma=True)
    if reciprocal_density:
        v.config_dict["KPOINTS"].update(
            {"reciprocal_density":reciprocal_density})
        v = DictSet(structure, v.config_dict)
    fws = []

    fws.append(OptimizeFW(structure=structure,
                          vasp_input_set=v,
                          vasp_cmd=vasp_cmd,
                          db_file=db_file))

    deformations = []
    # Generate deformations
    for ind in [(0, 0), (1, 1), (2, 2)]:
        for amount in norm_deformations:
            defo = Deformation.from_index_amount(ind, amount)
            deformations.append(defo)

    for ind in [(0, 1), (0, 2), (1, 2)]:
        for amount in shear_deformations:
            defo = Deformation.from_index_amount(ind, amount)
            deformations.append(defo)

    def_vasp_params = {"user_incar_settings": v.incar.as_dict()}
    def_vasp_params["user_incar_settings"].update({"ISIF":2,"ISTART":1})

    if reciprocal_density:
        def_vasp_params.update(
            {"reciprocal_density":reciprocal_density})
    
    for deformation in deformations:
        # TODO: Maybe should be more general, needing to specify
        #   the vasp input set with a string is a bit unwieldy
        #   for complete customization of the INCAR parameters
        fw = TransmuterFW(name="elastic deformation",
                          structure=structure,
                          transformations=['DeformStructureTransformation'],
                          transformation_params=[
                              {"deformation": deformation.tolist()}],
                          copy_vasp_outputs=True,
                          db_file=db_file,
                          vasp_cmd=vasp_cmd,
                          parents=fws[0],
                          vasp_input_params = def_vasp_params
                         )
        fw.spec['_tasks'].append(
            PassStressStrainData(deformation=deformation.tolist()).to_dict())
        fws.append(fw)
    
    fws.append(Firework(AnalyzeStressStrainData(structure=structure, 
                                                db_file=db_file),
                        name="Analyze Elastic Data", parents=fws[1:],
                        spec = {"_allow_fizzled_parents":True}))

    wfname = "{}:{}".format(structure.composition.reduced_formula,
                            "elastic constants")
    return Workflow(fws, name=wfname)

if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest

    structure = PymatgenTest.get_structure("Si")
    wf = get_wf_elastic_constant(structure)
