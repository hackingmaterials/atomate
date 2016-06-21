# coding: utf-8

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
from decimal import Decimal

import numpy as np
from fireworks import FireTaskBase, Firework, FWAction, Workflow
from fireworks.features.dupefinder import DupeFinderBase
from fireworks.utilities.fw_utilities import explicit_serialize

from matmethods.vasp.fireworks.core import OptimizeFW, TransmuterFW
from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.analysis.elasticity.strain import Deformation, IndependentStrain
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.sets import MPVaspInputSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.standard_transformations import \
    DeformStructureTransformation

"""
This module defines the elastic workflow
"""


__author__ = 'Shyam Dwaraknath, Joseph Montoya'
__email__ = 'shyamd@lbl.gov, montoyjh@lbl.gov'


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
#                            'fw_id': fw_spec['fw_id']}
        deformations.append(deformation_dict)
        return FWAction(mod_spec=[{'_push_all': {'deformations': deformations}}])


@explicit_serialize
class AnalyzeStressStrainData(FireTaskBase):
    """
    Passes the stress and deformation for an elastic deformation calculation

    Required params:
        deformation: the deformation gradient used in the elastic analysis.
    """

    required_params = ['structure']
    optional_params = ['db_file']

    def run_task(self, fw_spec):
        deformations = fw_spec['deformations']
        d = {"analysis": {}, "error": [], "warning": [],
             'deformation_tasks': {}, 'structure': self['structure']}
        stress_dict = {}

        dtypes = set()
        for deformation in deformations:
            defo = deformation['deformation']
            d_ind = np.nonzero(defo - np.eye(3))
            delta = Decimal((defo - np.eye(3))[d_ind][0])
            # Normal deformation
            if d_ind[0] == d_ind[1]:
                dtype = "_".join(["d", str(d_ind[0][0]),
                                  "{:.0e}".format(delta)])
                dtypes.add("_".join(["d", str(d_ind[0][0])]))
            # Shear deformation
            else:
                dtype = "_".join(["s", str(d_ind[0] + d_ind[1]),
                                  "{:.0e}".format(delta)])
                dtypes.add("_".join(["s", str(d_ind[0] + d_ind[1])]))

            strain = IndependentStrain(defo)
            stress = Stress(deformation['stress'])
            d["deformation_tasks"][dtype] = {'deformation_matrix': defo,
                                             'strain': strain.tolist(),
#                                             'task_id': deformation['task_id'],
                                             'stress': deformation['stress']}
            stress_dict[strain] = stress

        # DETERMINE IF WE HAVE 6 "UNIQUE" deformations
        if len(dtypes) == 6:
            # Perform Elastic tensor fitting and analysis
            result = ElasticTensor.from_stress_dict(stress_dict)
            d["elastic_tensor"] = result.tolist()
            kg_average = result.kg_average
            d.update({"K_Voigt": kg_average[0], "G_Voigt": kg_average[1],
                      "K_Reuss": kg_average[2], "G_Reuss": kg_average[3],
                      "K_Voigt_Reuss_Hill": kg_average[4],
                      "G_Voigt_Reuss_Hill": kg_average[5]})
            d["universal_anisotropy"] = result.universal_anisotropy
            d["homogeneous_poisson"] = result.homogeneous_poisson

            # Perform filter checks
            symm_t = result.symmetrized
            d["symmetrized_tensor"] = symm_t.tolist()
            d["analysis"]["not_rare_earth"] = True
            for s in self['structure'].species:
                if s.is_rare_earth_metal:
                    d["analysis"]["not_rare_earth"] = False
            eigvals = np.linalg.eigvals(symm_t)
            eig_positive = np.all((eigvals > 0) & np.isreal(eigvals))
            d["analysis"]["eigval_positive"] = bool(eig_positive)
            c11 = symm_t[0][0]
            c12 = symm_t[0][1]
            c13 = symm_t[0][2]
            c23 = symm_t[1][2]

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

            # TODO:
            # JHM: eventually we can reintroduce the IEEE conversion
            #       but as of now it's not being used, and it should
            #       be in pymatgen
            """
            # IEEE Conversion
            try:
                ieee_tensor = IEEE_conversion.get_ieee_tensor(struct_final, result)
                d["elastic_tensor_IEEE"] = ieee_tensor[0].tolist()
                d["analysis"]["IEEE"] = True
            except Exception as e:
                d["elastic_tensor_IEEE"] = None
                d["analysis"]["IEEE"] = False
                d["error"].append("Unable to get IEEE tensor: {}".format(e))
            """
            # Add thermal properties
            nsites = self['struct'].num_sites
            volume = self['struct'].volume
            natoms = self['struct'].composition.num_atoms
            weight = self['struct'].composition.weight
            num_density = 1e30 * nsites / volume
            mass_density = 1.6605e3 * nsites * volume * weight / \
                (natoms * volume)
            tot_mass = sum([e.atomic_mass for e in calc_struct.species])
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
                            "youngs_modulus": y_mod,
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
            d['state'] = "Fewer than 6 unique deformations"
            return FWAction()

        if d["error"]:
            raise ValueError("Elastic analysis failed: {}".format(d["error"]))
        elif d["analysis"]["filter_pass"]:
            d["state"] = "successful"
        else:
            d["state"] = "filter_failed"

        if not db_file:
            with open("elasticity.json", "w") as f:
                f.write(json.dumps(d, default=DATETIME_HANDLER))
        else:
            db_config = get_settings(db_file)
            db = MMDb(host=db_config["host"], port=db_config["port"],
                      database=db_config["database"],
                      user=db_config.get("admin_user"),
                      password=db_config.get("admin_password"),
                      collection='elasticity')
            db.collection.insert(d)
            logger.info("Finished analyzing elasticity data")


class DupeFinderEquivelantDeformation(DupeFinderBase):
    """
    This DupeFinder finds symmetrically equiavelant deformations
    """

    _fw_name = 'DupeFinderEquivelantDeformation'

    def verify(self, spec1, spec2):
        # Same type of firework and composition
        if spec1['name'] is not spec2['name']:
            return False

        # Find the write io task and check structure
        writeio1 = [task for task in spec1['_task'] if task['_fw_name']
                    is "{{matmethods.vasp.firetasks.write_inputs.WriteTransmutedStructureIOSet}}"]
        writeio2 = [task for task in spec2['_task'] if task['_fw_name']
                    is "{{matmethods.vasp.firetasks.write_inputs.WriteTransmutedStructureIOSet}}"]

        if writeio1['structure'] is not writeio2['structure']:
            return False

        # Are the deformations symmetrically equivelant?
        sga = SpacegroupAnalyzer(writeio1['structure'], tol=0.1)
        symm_ops = sga.get_symmetry_operations(cartesian=True)
        if not equivelant_by_symmetry(writeio1['transformation_params']['deformation'],
                                      writeio2['transformation_params']['deformation'],
                                      symm_ops):
            return False

        return True

    def query(self, spec):
        # QUERY FOR A EQUIVELANT TRANSFORM STRUCTURE FIREWORK
        # TODO: How do we provide a replacement that is the transformed stress/strain to match

        # Query for task name and structure
        return {'name': spec['name']}


def get_wf_elastic_constant(structure, norm_deformations=[-0.01, -0.005, 0.005, 0.01],
                            shear_deformations=[-0.08, -0.04, 0.04, 0.08], vasp_input_set=None,
                            vasp_cmd="vasp", db_file=None):
    """
    Returns a workflow to calculate elastic consants. This workflow is dynamic

    Firework 1 : write vasp input set for structural relaxation,
                 run vasp,
                 pass run location,
                 database insertion.

    Firework 2: Determine which deformations to run based on symmetry if needed
                TODO: How to turn off children selectively?

    Firework 3 - N: Optimize Deformed Structure

    Firework N: Analyze Deformations


    Args:
        structure (Structure): input structure to be optimized and run
        vasp_input_set (DictVaspInputSet): vasp input set.
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.

    Returns:
        Workflow
    """

    v = vasp_input_set or MPVaspInputSet(force_gamma=True)

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

    for deformation in deformations:
        fw = TransmuterFW(name="elastic deformation calculation",
                          structure=structure,
                          transformations=["DeformStructureTransformation"],
                          transformation_params=[
                              {"deformation": deformation.tolist()}],
                          copy_vasp_outputs=True,
                          db_file=db_file,
                          vasp_cmd=vasp_cmd,
                          parents=fws[0],
                          vasp_input_params={"reciprocal_density": 400,
                                             "user_incar_settings": {"ISIF": 2, "ENCUT": 700, "EDIFF": 0.000001}})

        fw.spec['_tasks'].append(PassStressStrainData(deformation=deformation.tolist()).to_dict())
        fws.append(fw)

    fws.append(Firework(AnalyzeStressStrainData(structure=structure),
                        name="Analyze Elastic Data", parents=fws[1:]))

    wfname = "{}:{}".format(structure.composition.reduced_formula, "elastic constants")
    return Workflow(fws, name=wfname)


def equivelant_by_symmetry(tensor1, tensor2, symm_ops, tolerance=1e-2):

    for sym in symm_ops:
        if (np.abs(symm_op.transform_tensor(tensor1) - tensor2) < tol).all():
            return True

    return False

if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest

    structure = PymatgenTest.get_structure("Si")
    wf = get_wf_elastic_constant(structure)
