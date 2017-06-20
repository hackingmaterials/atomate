# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import json
import os
from collections import defaultdict
from datetime import datetime

import numpy as np

from monty.json import MontyEncoder

from fireworks import FiretaskBase, FWAction, explicit_serialize
from fireworks.utilities.fw_serializers import DATETIME_HANDLER

from pymatgen import Structure
from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.analysis.elasticity.strain import IndependentStrain, Strain
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.electronic_structure.boltztrap import BoltztrapAnalyzer
from pymatgen.io.vasp.sets import get_vasprun_outcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from atomate.common.firetasks.glue_tasks import get_calc_loc
from atomate.utils.utils import env_chk, get_meta_from_structure
from atomate.utils.utils import get_logger
from atomate.vasp.database import VaspCalcDb
from atomate.vasp.drones import VaspDrone

__author__ = 'Anubhav Jain, Kiran Mathew, Shyam Dwaraknath'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov, shyamd@lbl.gov'

logger = get_logger(__name__)


@explicit_serialize
class VaspToDb(FiretaskBase):
    """
    Enter a VASP run into the database. Uses current directory unless you
    specify calc_dir or calc_loc.

    Optional params:
        calc_dir (str): path to dir (on current filesystem) that contains VASP
            output files. Default: use current working directory.
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
            search for the most recent calc_loc with the matching name
        parse_dos (bool): whether to parse the DOS and store in GridFS.
            Defaults to False.
        bandstructure_mode (str): Set to "uniform" for uniform band structure.
            Set to "line" for line mode. If not set, band structure will not
            be parsed.
        additional_fields (dict): dict of additional fields to add
        db_file (str): path to file containing the database credentials.
            Supports env_chk. Default: write data to JSON file.
        fw_spec_field (str): if set, will update the task doc with the contents
            of this key in the fw_spec.
        defuse_unsuccessful (bool): Defuses children fireworks if VASP run state
            is not "successful"; i.e. both electronic and ionic convergence are reached.
            Defaults to True.
    """
    optional_params = ["calc_dir", "calc_loc", "parse_dos", "bandstructure_mode",
                       "additional_fields", "db_file", "fw_spec_field", "defuse_unsuccessful"]

    def run_task(self, fw_spec):
        # get the directory that contains the VASP dir to parse
        calc_dir = os.getcwd()
        if "calc_dir" in self:
            calc_dir = self["calc_dir"]
        elif self.get("calc_loc"):
            calc_dir = get_calc_loc(self["calc_loc"], fw_spec["calc_locs"])["path"]

        # parse the VASP directory
        logger.info("PARSING DIRECTORY: {}".format(calc_dir))

        drone = VaspDrone(additional_fields=self.get("additional_fields"),
                          parse_dos=self.get("parse_dos", False), compress_dos=1,
                          bandstructure_mode=self.get("bandstructure_mode", False), compress_bs=1)

        # assimilate (i.e., parse)
        task_doc = drone.assimilate(calc_dir)

        # Check for additional keys to set based on the fw_spec
        if self.get("fw_spec_field"):
            task_doc.update(fw_spec[self.get("fw_spec_field")])

        # get the database connection
        db_file = env_chk(self.get('db_file'), fw_spec)

        # db insertion or taskdoc dump
        if not db_file:
            with open("task.json", "w") as f:
                f.write(json.dumps(task_doc, default=DATETIME_HANDLER))
        else:
            mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
            t_id = mmdb.insert_task(task_doc,
                                    parse_dos=self.get("parse_dos", False),
                                    parse_bs=bool(self.get("bandstructure_mode", False)))
            logger.info("Finished parsing with task_id: {}".format(t_id))

        if self.get("defuse_unsuccessful", True):
            defuse_children = (task_doc["state"] != "successful")
        else:
            defuse_children = False

        return FWAction(stored_data={"task_id": task_doc.get("task_id", None)},
                        defuse_children=defuse_children)


@explicit_serialize
class JsonToDb(FiretaskBase):
    """
    Insert the a JSON file (default: task.json) directly into the tasks database.
    Note that if the JSON file contains a "task_id" key, that task_id must not already be present
    in the tasks collection.

    Optional params:
        json_filename (str): name of the JSON file to insert (default: "task.json")
        db_file (str): path to file containing the database credentials. Supports env_chk.
        calc_dir (str): path to dir (on current filesystem) that contains VASP output files.
            Default: use current working directory.
    """
    optional_params = ["json_filename", "db_file", "calc_dir"]

    def run_task(self, fw_spec):

        ref_file = self.get("json_filename", "task.json")
        calc_dir = self.get("calc_dir", os.getcwd())
        with open(os.path.join(calc_dir, ref_file), "r") as fp:
            task_doc = json.load(fp)

        db_file = env_chk(self.get('db_file'), fw_spec)
        if not db_file:
            with open("task.json", "w") as f:
                f.write(json.dumps(task_doc, default=DATETIME_HANDLER))
        else:
            mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
            mmdb.insert(task_doc)


@explicit_serialize
class BoltztrapToDb(FiretaskBase):
    """
    Enter a BoltzTraP run into the database. Note that this assumes you are in a current dir
    that has the uniform band structure data with a sub-directory called "boltztrap" containing
    the BoltzTraP information.

    Optional params:
        db_file (str): path to file containing the database credentials.
            Supports env_chk. Default: write data to JSON file.
        hall_doping (bool): set True to retain hall_doping in dict
        additional_fields (dict): fields added to the document such as user-defined tags or name, ids, etc
    """

    optional_params = ["db_file", "hall_doping", "additional_fields"]

    def run_task(self, fw_spec):
        additional_fields = self.get("additional_fields", {})

        # pass the additional_fields first to avoid overriding BoltztrapAnalyzer items
        d = additional_fields.copy()

        btrap_dir = os.path.join(os.getcwd(), "boltztrap")
        d["boltztrap_dir"] = btrap_dir

        bta = BoltztrapAnalyzer.from_files(btrap_dir)
        d.update(bta.as_dict())
        d["scissor"] = bta.intrans["scissor"]

        # trim the output
        for x in ['cond', 'seebeck', 'kappa', 'hall', 'mu_steps', 'mu_doping', 'carrier_conc']:
            del d[x]

        if not self.get("hall_doping"):
            del d["hall_doping"]

        bandstructure_dir = os.getcwd()
        d["bandstructure_dir"] = bandstructure_dir

        # add the structure
        v, o = get_vasprun_outcar(bandstructure_dir, parse_eigen=False, parse_dos=False)
        structure = v.final_structure
        d["structure"] = structure.as_dict()
        d.update(get_meta_from_structure(structure))

        # add the spacegroup
        sg = SpacegroupAnalyzer(Structure.from_dict(d["structure"]), 0.1)
        d["spacegroup"] = {"symbol": sg.get_space_group_symbol(),
                           "number": sg.get_space_group_number(),
                           "point_group": sg.get_point_group_symbol(),
                           "source": "spglib",
                           "crystal_system": sg.get_crystal_system(),
                           "hall": sg.get_hall()}

        d["created_at"] = datetime.utcnow()

        db_file = env_chk(self.get('db_file'), fw_spec)

        if not db_file:
            del d["dos"]
            with open(os.path.join(btrap_dir, "boltztrap.json"), "w") as f:
                f.write(json.dumps(d, default=DATETIME_HANDLER))
        else:
            mmdb = VaspCalcDb.from_db_file(db_file, admin=True)

            # dos gets inserted into GridFS
            dos = json.dumps(d["dos"], cls=MontyEncoder)
            fsid, compression = mmdb.insert_gridfs(dos, collection="dos_boltztrap_fs",
                                                   compress=True)
            d["dos_boltztrap_fs_id"] = fsid
            del d["dos"]

            mmdb.db.boltztrap.insert(d)


@explicit_serialize
class ElasticTensorToDb(FiretaskBase):
    """
    Analyzes the stress/strain data of an elastic workflow to produce
    an elastic tensor and various other quantities.
    """

    required_params = ['structure']
    optional_params = ['db_file']

    def run_task(self, fw_spec):
        d = {"analysis": {},
             "deformation_tasks": fw_spec["deformation_tasks"],
             "initial_structure": self['structure'].as_dict()}

        # Get optimized structure
        calc_locs_opt = [cl for cl in fw_spec['calc_locs'] if 'optimize' in cl['name']]
        if calc_locs_opt:
            optimize_loc = calc_locs_opt[-1]['path']
            logger.info("Parsing initial optimization directory: {}".format(optimize_loc))
            drone = VaspDrone()
            optimize_doc = drone.assimilate(optimize_loc)
            opt_struct = Structure.from_dict(optimize_doc["calcs_reversed"][0]["output"]["structure"])
            d.update({"optimized_structure": opt_struct.as_dict()})

        # TODO: @montoyjh: does the below have anything to do with elastic tensor? If not, try
        # the more general fw_spec_field approach in the VaspToDb rather than hard-coding the
        # tags insertion here. -computron
        if fw_spec.get("tags", None):
            d["tags"] = fw_spec["tags"]

        results = fw_spec["deformation_tasks"].values()
        defos = [r["deformation_matrix"] for r in results]
        stresses = [r["stress"] for r in results]
        strains = np.array([Strain(r["strain"]).voigt for r in results])
        stress_dict = {IndependentStrain(defo) : Stress(stress) for defo, stress in zip(defos, stresses)}

        logger.info("Analyzing stress/strain data")
        # Determine if we have 6 unique deformations
        # TODO: @montoyjh: what if it's a cubic system? don't need 6. -computron
        if np.linalg.matrix_rank(strains) == 6:
            # Perform Elastic tensor fitting and analysis
            result = ElasticTensor.from_stress_dict(stress_dict)
            d["elastic_tensor"] = result.voigt.tolist()
            d.update(result.property_dict)

        else:
            raise ValueError("Fewer than 6 unique deformations")

        d["state"] = "successful"

        # Save analysis results in json or db
        db_file = env_chk(self.get('db_file'), fw_spec)
        if not db_file:
            with open("elasticity.json", "w") as f:
                f.write(json.dumps(d, default=DATETIME_HANDLER))
        else:
            db = VaspCalcDb.from_db_file(db_file, admin=True)
            db.collection = db.db["elasticity"]
            db.collection.insert_one(d)
            logger.info("Elastic analysis complete.")
        return FWAction()


@explicit_serialize
class RamanTensorToDb(FiretaskBase):
    """
    Raman susceptibility tensor for each mode = Finite difference derivative of the dielectric
        tensor wrt the displacement along that mode.
    See: 10.1103/PhysRevB.73.104304.
    The frequencies are in the units of cm^-1. To convert the frequency to THz: multiply by 0.1884.


    optional_params:
        db_file (str): path to the db file
    """

    optional_params = ["db_file"]

    def run_task(self, fw_spec):
        nm_eigenvecs = np.array(fw_spec["normalmodes"]["eigenvecs"])
        nm_eigenvals = np.array(fw_spec["normalmodes"]["eigenvals"])
        nm_norms = np.linalg.norm(nm_eigenvecs, axis=2)
        structure = fw_spec["normalmodes"]["structure"]
        masses = np.array([site.specie.data['Atomic mass'] for site in structure])
        nm_norms = nm_norms / np.sqrt(masses)  # eigenvectors in vasprun.xml are not divided by sqrt(M_i)
        # To get the actual eigenvals, the values read from vasprun.xml must be multiplied by -1.
        # frequency_i = sqrt(-e_i)
        # To convert the frequency to THZ: multiply sqrt(-e_i) by 15.633
        # To convert the frequency to cm^-1: multiply sqrt(-e_i) by 82.995
        nm_frequencies = np.sqrt(np.abs(nm_eigenvals)) * 82.995  # cm^-1

        d = {"structure": structure.as_dict(),
             "normalmodes": {"eigenvals": fw_spec["normalmodes"]["eigenvals"],
                             "eigenvecs": fw_spec["normalmodes"]["eigenvecs"]
                             },
             "frequencies": nm_frequencies.tolist()}

        # store the displacement & epsilon for each mode in a dictionary
        mode_disps = fw_spec["raman_epsilon"].keys()
        modes_eps_dict = defaultdict(list)
        for md in mode_disps:
            modes_eps_dict[fw_spec["raman_epsilon"][md]["mode"]].append(
                [fw_spec["raman_epsilon"][md]["displacement"],
                 fw_spec["raman_epsilon"][md]["epsilon"]])

        # raman tensor = finite difference derivative of epsilon wrt displacement.
        raman_tensor_dict = {}
        scale = np.sqrt(structure.volume/2.0) / 4.0 / np.pi
        for k, v in modes_eps_dict.items():
            raman_tensor = (np.array(v[0][1]) - np.array(v[1][1])) / (v[0][0] - v[1][0])
            # frequency in cm^-1
            omega = nm_frequencies[k]
            if nm_eigenvals[k] > 0:
                logger.warn("Mode: {} is UNSTABLE. Freq(cm^-1) = {}".format(k, -omega))
            raman_tensor = scale * raman_tensor * np.sum(nm_norms[k]) / np.sqrt(omega)
            raman_tensor_dict[str(k)] = raman_tensor.tolist()

        d["raman_tensor"] = raman_tensor_dict
        d["state"] = "successful"

        # store the results
        db_file = env_chk(self.get("db_file"), fw_spec)
        if not db_file:
            with open("raman.json", "w") as f:
                f.write(json.dumps(d, default=DATETIME_HANDLER))
        else:
            db = VaspCalcDb.from_db_file(db_file, admin=True)
            db.collection = db.db["raman"]
            db.collection.insert_one(d)
            logger.info("Raman tensor calculation complete.")
        return FWAction()


# TODO: @computron: this requires a "tasks" collection to proceed. Merits of changing to FW passing
# method? -computron
# TODO: @computron: even if you use the db-centric method, embed information in tags rather than
# task_label? This workflow likely requires review with its authors. -computron
@explicit_serialize
class GibbsAnalysisToDb(FiretaskBase):
    """
    Compute the quasi-harmonic gibbs free energy. There are 2 options available for the
    quasi-harmonic approximation (set via 'qha_type' parameter):
    1. use the phonopy package quasi-harmonic approximation interface or
    2. use the debye model.
    Note: Instead of relying on fw_spec, this task gets the required data directly from the
    tasks collection for processing. The summary dict is written to 'gibbs.json' file.

    required_params:
        tag (str): unique tag appended to the task labels in other fireworks so that all the
            required data can be queried directly from the database.
        db_file (str): path to the db file

    optional_params:
        qha_type(str): quasi-harmonic approximation type: "debye_model" or "phonopy",
            default is "debye_model"
        t_min (float): min temperature
        t_step (float): temperature step
        t_max (float): max temperature
        mesh (list/tuple): reciprocal space density
        eos (str): equation of state used for fitting the energies and the volumes.
            options supported by phonopy: "vinet", "murnaghan", "birch_murnaghan".
        pressure (float): in GPa, optional.
        metadata (dict): meta data

    """

    required_params = ["tag", "db_file"]
    optional_params = ["qha_type", "t_min", "t_step", "t_max", "mesh", "eos", "pressure", "poisson",
                       "metadata"]

    def run_task(self, fw_spec):

        gibbs_dict = {}

        tag = self["tag"]
        t_step = self.get("t_step", 10)
        t_min = self.get("t_min", 0)
        t_max = self.get("t_max", 1000)
        mesh = self.get("mesh", [20, 20, 20])
        eos = self.get("eos", "vinet")
        qha_type = self.get("qha_type", "debye_model")
        pressure = self.get("pressure", 0.0)
        poisson = self.get("poisson", 0.25)
        gibbs_dict["metadata"] = self.get("metadata", {})


        db_file = env_chk(self.get("db_file"), fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        # get the optimized structure
        d = mmdb.collection.find_one({"task_label": "{} structure optimization".format(tag)},
                                     {"calcs_reversed": 1})
        structure = Structure.from_dict(d["calcs_reversed"][-1]["output"]['structure'])
        gibbs_dict["structure"] = structure.as_dict()

        # get the data(energy, volume, force constant) from the deformation runs
        docs = mmdb.collection.find({"task_label": {"$regex": "{} gibbs*".format(tag)},
                                     "formula_pretty": structure.composition.reduced_formula},
                                    {"calcs_reversed": 1})
        energies = []
        volumes = []
        force_constants = []
        for d in docs:
            s = Structure.from_dict(d["calcs_reversed"][-1]["output"]['structure'])
            energies.append(d["calcs_reversed"][-1]["output"]['energy'])
            if qha_type not in ["debye_model"]:
                force_constants.append(d["calcs_reversed"][-1]["output"]['force_constants'])
            volumes.append(s.volume)
        gibbs_dict["energies"] = energies
        gibbs_dict["volumes"] = volumes
        if qha_type not in ["debye_model"]:
            gibbs_dict["force_constants"] = force_constants

        try:
            # use quasi-harmonic debye approximation
            if qha_type in ["debye_model"]:

                from pymatgen.analysis.quasiharmonic import QuasiharmonicDebyeApprox

                qhda = QuasiharmonicDebyeApprox(energies, volumes, structure, t_min, t_step, t_max,
                                                eos, pressure=pressure, poisson=poisson)
                gibbs_dict.update(qhda.get_summary_dict())
                gibbs_dict["success"] = True

            # use the phonopy interface
            else:

                from atomate.vasp.analysis.phonopy import get_phonopy_gibbs

                G, T = get_phonopy_gibbs(energies, volumes, force_constants, structure, t_min,
                                         t_step, t_max, mesh, eos, pressure)
                gibbs_dict["gibbs_free_energy"] = G
                gibbs_dict["temperatures"] = T
                gibbs_dict["success"] = True

        # quasi-harmonic analysis failed, set the flag to false
        except:
            import traceback

            logger.warn("Quasi-harmonic analysis failed!")
            gibbs_dict["success"] = False
            gibbs_dict["traceback"] = traceback.format_exc()
            metadata.update({"task_label_tag": tag})
            gibbs_dict["metadata"] = metadata
            gibbs_dict["created_at"] = datetime.utcnow()

        # TODO: @matk86: add a list of task_ids that were used to construct the analysis to DB?
        # -computron
        if not db_file:
            dump_file = "gibbs.json"
            logger.info("Dumping the analysis summary to {}".format(dump_file))
            with open(dump_file, "w") as f:
                f.write(json.dumps(gibbs_dict, default=DATETIME_HANDLER))
        else:
            coll = mmdb.db["gibbs_tasks"]
            coll.insert_one(gibbs_dict)

        logger.info("Gibbs free energy calculation complete.")

        if not gibbs_dict["success"]:
            return FWAction(defuse_children=True)


# TODO: @computron: review method of data passing with the workflow authors. -computron
@explicit_serialize
class FitEOSToDb(FiretaskBase):
    """
    Retrieve the energy and volume data and fit it to the given equation of state. The summary dict
    is written to 'bulk_modulus.json' file.

    Required parameters:
        tag (str): unique tag appended to the task labels in other fireworks so that all the
            required data can be queried directly from the database.
        db_file (str): path to the db file

    Optional parameters:
        to_db (bool): if True, the data will be inserted to "eos" collection; otherwise, dumped to a .json file.
        eos (str): equation of state used for fitting the energies and the volumes.
            options supported by pymatgen: "quadratic", "murnaghan", "birch", "birch_murnaghan",
            "pourier_tarantola", "vinet", "deltafactor". Default: "vinet"
    """

    required_params = ["tag", "db_file"]
    optional_params = ["to_db", "eos"]

    def run_task(self, fw_spec):

        from pymatgen.analysis.eos import EOS

        eos = self.get("eos", "vinet")
        tag = self["tag"]
        db_file = env_chk(self.get("db_file"), fw_spec)
        summary_dict = {"eos": eos}
        to_db = self.get("to_db", True)

        # collect and store task_id of all related tasks to make unique links with "tasks" collection
        all_task_ids = []

        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        # get the optimized structure
        d = mmdb.collection.find_one({"task_label": "{} structure optimization".format(tag)})
        all_task_ids.append(d["task_id"])
        structure = Structure.from_dict(d["calcs_reversed"][-1]["output"]['structure'])
        summary_dict["structure"] = structure.as_dict()

        # get the data(energy, volume, force constant) from the deformation runs
        docs = mmdb.collection.find({"task_label": {"$regex": "{} bulk_modulus*".format(tag)},
                                     "formula_pretty": structure.composition.reduced_formula})
        energies = []
        volumes = []
        for d in docs:
            s = Structure.from_dict(d["calcs_reversed"][-1]["output"]['structure'])
            energies.append(d["calcs_reversed"][-1]["output"]['energy'])
            volumes.append(s.volume)
            all_task_ids.append(d["task_id"])
        summary_dict["energies"] = energies
        summary_dict["volumes"] = volumes
        summary_dict["all_task_ids"] = all_task_ids

        # fit the equation of state
        eos = EOS(eos)
        eos_fit = eos.fit(volumes, energies)
        summary_dict["bulk_modulus"] = eos_fit.b0_GPa

        # TODO: find a better way for passing tags of the entire workflow to db - albalu
        if fw_spec.get("tags", None):
            summary_dict["tags"] = fw_spec["tags"]
        summary_dict["results"] = dict(eos_fit.results)
        summary_dict["created_at"] = datetime.utcnow()

        # db_file itself is required but the user can choose to pass the results to db or not
        if to_db:
            mmdb.collection = mmdb.db["eos"]
            mmdb.collection.insert_one(summary_dict)
        else:
            with open("bulk_modulus.json", "w") as f:
                f.write(json.dumps(summary_dict, default=DATETIME_HANDLER))

        # TODO: @matk86 - there needs to be a builder to put it into materials collection... -computron
        logger.info("Bulk modulus calculation complete.")


# TODO: @computron: review method of data passing with the workflow authors. -computron
# TODO: insert to db. -matk
@explicit_serialize
class ThermalExpansionCoeffToDb(FiretaskBase):
    """
    Compute the quasi-harmonic thermal expansion coefficient using phonopy.

    required_params:
        tag (str): unique tag appended to the task labels in other fireworks so that all the
            required data can be queried directly from the database.
        db_file (str): path to the db file

    optional_params:
        t_min (float): min temperature
        t_step (float): temperature step
        t_max (float): max temperature
        mesh (list/tuple): reciprocal space density
        eos (str): equation of state used for fitting the energies and the volumes.
            options supported by phonopy: "vinet" (default), "murnaghan", "birch_murnaghan".
        pressure (float): in GPa, optional.
    """

    required_params = ["tag", "db_file"]
    optional_params = ["t_min", "t_step", "t_max", "mesh", "eos", "pressure"]

    def run_task(self, fw_spec):

        from atomate.vasp.analysis.phonopy import get_phonopy_thermal_expansion

        tag = self["tag"]
        db_file = env_chk(self.get("db_file"), fw_spec)
        t_step = self.get("t_step", 10)
        t_min = self.get("t_min", 0)
        t_max = self.get("t_max", 1000)
        mesh = self.get("mesh", [20, 20, 20])
        eos = self.get("eos", "vinet")
        pressure = self.get("pressure", 0.0)
        summary_dict = {}

        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        # get the optimized structure
        d = mmdb.collection.find_one({"task_label": "{} structure optimization".format(tag)})
        structure = Structure.from_dict(d["calcs_reversed"][-1]["output"]['structure'])
        summary_dict["structure"] = structure.as_dict()

        # get the data(energy, volume, force constant) from the deformation runs
        docs = mmdb.collection.find({"task_label": {"$regex": "{} thermal_expansion*".format(tag)},
                                     "formula_pretty": structure.composition.reduced_formula})
        energies = []
        volumes = []
        force_constants = []
        for d in docs:
            s = Structure.from_dict(d["calcs_reversed"][-1]["output"]['structure'])
            energies.append(d["calcs_reversed"][-1]["output"]['energy'])
            volumes.append(s.volume)
            force_constants.append(d["calcs_reversed"][-1]["output"]['force_constants'])
        summary_dict["energies"] = energies
        summary_dict["volumes"] = volumes
        summary_dict["force_constants"] = force_constants

        alpha, T = get_phonopy_thermal_expansion(energies, volumes, force_constants, structure,
                                                 t_min, t_step, t_max, mesh, eos, pressure)

        summary_dict["alpha"] = alpha
        summary_dict["T"] = T

        with open("thermal_expansion.json", "w") as f:
            f.write(json.dumps(summary_dict, default=DATETIME_HANDLER))

        # TODO: @matk86 - there needs to be a way to insert this into a database! And also
        # a builder to put it into materials collection... -computron
        logger.info("Thermal expansion coefficient calculation complete.")


# the following definitions for backward compatibility
class VaspToDbTask(VaspToDb):
    pass


class JsonToDbTask(JsonToDb):
    pass


class BoltztrapToDBTask(BoltztrapToDb):
    pass


class ElasticTensorToDbTask(ElasticTensorToDb):
    pass


class RamanSusceptibilityTensorToDbTask(RamanTensorToDb):
    pass


class GibbsFreeEnergyTask(GibbsAnalysisToDb):
    pass


class FitEquationOfStateTask(FitEOSToDb):
    pass


class ThermalExpansionCoeffTask(ThermalExpansionCoeffToDb):
    pass
