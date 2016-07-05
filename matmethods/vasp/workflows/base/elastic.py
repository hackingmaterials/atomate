# coding: utf-8

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import json
import os
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
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen import Structure
from matgendb.util import get_settings

from pymatgen.transformations.standard_transformations import \
    DeformStructureTransformation

"""
This module defines the elastic workflow
"""

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
        deformations = list(fw_spec.get("deformations", []))
        v = Vasprun('vasprun.xml.gz')
        stress = v.ionic_steps[-1]['stress']
        deformation_dict = {'deformation': self['deformation'],
                            'stress': stress}
        deformations.append(deformation_dict)
        return FWAction(mod_spec=[{'_push_all': {'deformations': deformations}}])


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
        optimize_loc = fw_spec["calc_locs"][0]["path"]
        logger.info("PARSING INITIAL OPTIMIZATION DIRECTORY: {}".format(optimize_loc))
        drone = VaspDrone()
        optimize_doc = drone.assimilate(optimize_loc)
        opt_struct = Structure.from_dict(
            optimize_doc["calcs_reversed"][0]["output"]["structure"])
        
        deformations = fw_spec['deformations']
        d = {"analysis": {}, "deformation_tasks": {},
             "initial_structure": self['structure'].as_dict(), 
             "optimized_structure": opt_struct.as_dict()}
        stress_dict = {}

        dtypes = []
        for deformation in deformations:
            defo = deformation['deformation']
            d_ind = np.nonzero(defo - np.eye(3))
            delta = Decimal((defo - np.eye(3))[d_ind][0])
            # Shorthand is d_X_V, X is voigt index, V is value
            dtype = "_".join(["d", str(reverse_voigt_map[d_ind][0]),
                              "{:.0e}".format(delta)])
            strain = IndependentStrain(defo)
            stress = Stress(deformation['stress'])
            d["deformation_tasks"][dtype] = {'deformation_matrix': defo,
                                             'strain': strain.tolist(),
                                             'stress': deformation['stress']}
            dtypes.append(dtype)
            stress_dict[strain] = stress

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

            # Perform filter checks
            symm_t = result.voigt_symmetrized
            fit_t = result.fit_to_structure(opt_struct)
            d["analysis"]["fit_to_structure"] = \
                    bool(result.is_fit_to_structure(opt_struct))
            d["symmetrized_tensor"] = symm_t.voigt.tolist()
            d["structure_fit_tensor"] = fit_t.voigt.tolist()
            d["analysis"]["not_rare_earth"] = True
            for s in opt_struct.species:
                if s.is_rare_earth_metal:
                    d["analysis"]["not_rare_earth"] = False
            eigvals = np.linalg.eigvals(symm_t.voigt)
            eig_positive = bool((eigvals > 0).all() 
                                and np.isreal(eigvals).all())
            d["analysis"]["eigval_positive"] = eig_positive
            c11 = symm_t.voigt[0][0]
            c12 = symm_t.voigt[0][1]
            c13 = symm_t.voigt[0][2]
            c23 = symm_t.voigt[1][2]

            d["analysis"]["c11_c12"] = not (abs((c11 - c12) / c11) < 0.05
                                            or c11 < c12)
            d["analysis"]["c11_c13"] = not (abs((c11 - c13) / c11) < 0.05
                                            or c11 < c13)
            d["analysis"]["c11_c23"] = not (abs((c11 - c23) / c11) < 0.1
                                            or c11 < c23)
            d["analysis"]["K_R"] = not (d["K_Reuss"] < 2)
            d["analysis"]["G_R"] = not (d["G_Reuss"] < 2)
            d["analysis"]["K_V"] = not (d["K_Voigt"] < 2)
            d["analysis"]["G_V"] = not (d["G_Voigt"] < 2)
            filter_state = np.all(d["analysis"].values())
            d["analysis"]["filter_pass"] = bool(filter_state)
            d["analysis"]["eigval"] = list(eigvals)
            
            # IEEE tensor
            ieee_tensor = result.convert_to_ieee(opt_struct)
            fit_ieee_tensor = fit_t.convert_to_ieee(opt_struct)
            d["elastic_tensor_IEEE"] = ieee_tensor.voigt.tolist()
            d["structure_fit_IEEE"] = fit_ieee_tensor.voigt.tolist()
            
            # Thermal properties
            nsites = opt_struct.num_sites
            volume = opt_struct.volume
            natoms = opt_struct.composition.num_atoms
            weight = opt_struct.composition.weight
            num_density = 1e30 * nsites / volume
            mass_density = 1.6605e3 * nsites * volume * weight / \
                (natoms * volume)
            tot_mass = sum([e.atomic_mass for e in opt_struct.species])
            avg_mass = 1.6605e-27 * tot_mass / natoms
            y_mod = 9e9 * result.k_vrh * result.g_vrh / \
                (3. * result.k_vrh * result.g_vrh)
            trans_v = 1e9 * result.k_vrh / mass_density**0.5
            long_v = 1e9 * result.k_vrh + \
                4. / 3. * result.g_vrh / mass_density**0.5
            clarke = 0.87 * 1.3806e-23 * avg_mass**(-2. / 3.) * \
                mass_density**(1. / 6.) * y_mod**0.5
            cahill = 1.3806e-23 / 2.48 * num_density**(2. / 3.) * long_v + \
                2 * trans_v
            snyder_ac = 0.38483 * avg_mass * \
                (long_v + 2. / 3. * trans_v)**3. / \
                (300. * num_density**(-2. / 3.) * nsites**(1. / 3.))
            snyder_opt = 1.66914e-23 * (long_v + 2. / 3. * trans_v) / \
                num_density**(-2. / 3.) * \
                (1 - nsites**(-1. / 3.))
            snyder_total = snyder_ac + snyder_opt
            debye = 2.489e-11 * avg_mass**(-1. / 3.) * \
                mass_density**(-1. / 6.) * y_mod**0.5

            d["thermal"] = {"num_density": num_density,
                            "mass_density": mass_density,
                            "avg_mass": avg_mass,
                            "num_atom_per_unit_formula": natoms,
                            "youngs": y_mod,
                            "trans_velocity": trans_v,
                            "long_velocity": long_v,
                            "clarke": clarke,
                            "cahill": cahill,
                            "snyder_acou_300K": snyder_ac,
                            "snyder_opt": snyder_opt,
                            "snyder_total": snyder_total,
                            "debye": debye
                            }
        else:
            raise ValueError("Fewer than 6 unique deformations")

        d["state"] = "successful"

        # get the database connection
        db_file = env_chk(self.get('db_file'), fw_spec)

        if not db_file:
            with open("elasticity.json", "w") as f:
                f.write(json.dumps(d, default=DATETIME_HANDLER))
        else:
            db = MMDb.from_db_file(db_file, admin=True)
            db.collection = db.db["elasticity"]
            import pdb; pdb.set_trace()
            db.collection.insert_one(d)
            logger.info("ELASTIC ANALYSIS COMPLETE")


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
    fws = []

    fws.append(OptimizeFW(structure=structure,
                          vasp_input_set=vasp_input_set,
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

    def_vasp_params = {"user_incar_settings":{"ISIF":2, "IBRION":2, 
                                              "NSW":99, "LAECHG":False,
                                              "LHVAR":False, "ALGO":"Fast",
                                              "LWAVE":False}}
    if reciprocal_density:
        deformation_vasp_params.update(
            {"reciprocal_density":reciprocal_density})
    
    for deformation in deformations:
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
