# coding: utf-8


import json
import os
import re
from collections import defaultdict
from datetime import datetime

import numpy as np

from monty.json import MontyEncoder, jsanitize
from pydash.objects import has, get

from atomate.vasp.config import DEFUSE_UNSUCCESSFUL
from fireworks import FiretaskBase, FWAction, explicit_serialize
from fireworks.utilities.fw_serializers import DATETIME_HANDLER

from pymatgen import Structure
from pymatgen.analysis.elasticity.elastic import ElasticTensor, ElasticTensorExpansion
from pymatgen.analysis.elasticity.strain import Strain, Deformation
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.electronic_structure.boltztrap import BoltztrapAnalyzer
from pymatgen.io.vasp.sets import get_vasprun_outcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.ferroelectricity.polarization import Polarization, get_total_ionic_dipole, \
    EnergyTrend
from pymatgen.analysis.magnetism import CollinearMagneticStructureAnalyzer, Ordering, magnetic_deformation
from pymatgen.command_line.bader_caller import bader_analysis_from_path

from atomate.common.firetasks.glue_tasks import get_calc_loc
from atomate.utils.utils import env_chk, get_meta_from_structure
from atomate.utils.utils import get_logger
from atomate.vasp.database import VaspCalcDb
from atomate.vasp.drones import VaspDrone
from atomate.vasp.config import STORE_VOLUMETRIC_DATA

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
        parse_potcar_file (bool): Whether to parse the potcar file. Defaults to
            True.
        bandstructure_mode (str): Set to "uniform" for uniform band structure.
            Set to "line" for line mode. If not set, band structure will not
            be parsed.
        additional_fields (dict): dict of additional fields to add
        db_file (str): path to file containing the database credentials.
            Supports env_chk. Default: write data to JSON file.
        fw_spec_field (str): if set, will update the task doc with the contents
            of this key in the fw_spec.
        defuse_unsuccessful (bool): this is a three-way toggle on what to do if
            your job looks OK, but is actually unconverged (either electronic or
            ionic). True -> mark job as COMPLETED, but defuse children.
            False --> do nothing, continue with workflow as normal. "fizzle"
            --> throw an error (mark this job as FIZZLED)
        task_fields_to_push (dict): if set, will update the next Firework/Firetask
            spec using fields from the task document.
            Format: {key : path} -> fw.spec[key] = task_doc[path]
            The path is a full mongo-style path so subdocuments can be referneced
            using dot notation and array keys can be referenced using the index.
            E.g "calcs_reversed.0.output.outar.run_stats"
    """
    optional_params = ["calc_dir", "calc_loc", "parse_dos", "bandstructure_mode",
                       "additional_fields", "db_file", "fw_spec_field", "defuse_unsuccessful",
                       "task_fields_to_push", "parse_chgcar", "parse_aeccar",
                       "parse_potcar_file",
                       "store_volumetric_data"]

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
                          parse_dos=self.get("parse_dos", False),
                          parse_potcar_file=self.get("parse_potcar_file", True),
                          bandstructure_mode=self.get("bandstructure_mode", False),
                          parse_chgcar=self.get("parse_chgcar", False),  # deprecated
                          parse_aeccar=self.get("parse_aeccar", False),  # deprecated
                          store_volumetric_data=self.get("store_volumetric_data", STORE_VOLUMETRIC_DATA))

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
            t_id = mmdb.insert_task(
                task_doc, use_gridfs=self.get("parse_dos", False)
                or bool(self.get("bandstructure_mode", False))
                or self.get("parse_chgcar", False)  # deprecated
                or self.get("parse_aeccar", False)  # deprecated
                or bool(self.get("store_volumetric_data", STORE_VOLUMETRIC_DATA)))
            logger.info("Finished parsing with task_id: {}".format(t_id))

        defuse_children = False
        if task_doc["state"] != "successful":
            defuse_unsuccessful = self.get("defuse_unsuccessful",
                                           DEFUSE_UNSUCCESSFUL)
            if defuse_unsuccessful is True:
                defuse_children = True
            elif defuse_unsuccessful is False:
                pass
            elif defuse_unsuccessful == "fizzle":
                raise RuntimeError(
                    "VaspToDb indicates that job is not successful "
                    "(perhaps your job did not converge within the "
                    "limit of electronic/ionic iterations)!")
            else:
                raise RuntimeError("Unknown option for defuse_unsuccessful: "
                                   "{}".format(defuse_unsuccessful))

        task_fields_to_push = self.get("task_fields_to_push", None)
        update_spec = {}
        if task_fields_to_push:
            if isinstance(task_fields_to_push, dict):
                for key, path_in_task_doc in task_fields_to_push.items():
                    if has(task_doc, path_in_task_doc):
                        update_spec[key] = get(task_doc, path_in_task_doc)
                    else:
                        logger.warning("Could not find {} in task document. Unable to push to next firetask/firework".format(path_in_task_doc))
            else:
                raise RuntimeError("Inappropriate type {} for task_fields_to_push. It must be a "
                                   "dictionary of format: {key: path} where key refers to a field "
                                   "in the spec and path is a full mongo-style path to a "
                                   "field in the task document".format(type(task_fields_to_push)))

        return FWAction(stored_data={"task_id": task_doc.get("task_id", None)},
                        defuse_children=defuse_children, update_spec=update_spec)


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
        d["formula_pretty"] = structure.composition.reduced_formula
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

    Required params:
        structure (Structure): structure to use for symmetrization,
            input structure.  If an optimization was used, will
            look for relaxed structure in calc locs

    Optional params:
        db_file (str): path to file containing the database credentials.
            Supports env_chk. Default: write data to JSON file.
        order (int): order of fit to perform
        fw_spec_field (str): if set, will update the task doc with the contents
            of this key in the fw_spec.
        fitting_method (str): if set, will use one of the specified
            fitting methods from pymatgen.  Supported methods are
            "independent", "pseudoinverse", and "finite_difference."
            Note that order 3 and higher required finite difference
            fitting, and will override.
    """

    required_params = ['structure']
    optional_params = ['db_file', 'order', 'fw_spec_field', 'fitting_method']

    def run_task(self, fw_spec):
        ref_struct = self['structure']
        d = {
            "analysis": {},
            "initial_structure": self['structure'].as_dict()
        }

        # Get optimized structure
        calc_locs_opt = [cl for cl in fw_spec.get('calc_locs', []) if 'optimiz' in cl['name']]
        if calc_locs_opt:
            optimize_loc = calc_locs_opt[-1]['path']
            logger.info("Parsing initial optimization directory: {}".format(optimize_loc))
            drone = VaspDrone()
            optimize_doc = drone.assimilate(optimize_loc)
            opt_struct = Structure.from_dict(optimize_doc["calcs_reversed"][0]["output"]["structure"])
            d.update({"optimized_structure": opt_struct.as_dict()})
            ref_struct = opt_struct
            eq_stress = -0.1*Stress(optimize_doc["calcs_reversed"][0]["output"]["ionic_steps"][-1]["stress"])
        else:
            eq_stress = None

        if self.get("fw_spec_field"):
            d.update({self.get("fw_spec_field"): fw_spec.get(self.get("fw_spec_field"))})

        # Get the stresses, strains, deformations from deformation tasks
        defo_dicts = fw_spec["deformation_tasks"].values()
        stresses, strains, deformations = [], [], []
        for defo_dict in defo_dicts:
            stresses.append(Stress(defo_dict["stress"]))
            strains.append(Strain(defo_dict["strain"]))
            deformations.append(Deformation(defo_dict["deformation_matrix"]))
            # Add derived stresses and strains if symmops is present
            for symmop in defo_dict.get("symmops", []):
                stresses.append(Stress(defo_dict["stress"]).transform(symmop))
                strains.append(Strain(defo_dict["strain"]).transform(symmop))
                deformations.append(Deformation(defo_dict["deformation_matrix"]).transform(symmop))

        stresses = [-0.1*s for s in stresses]
        pk_stresses = [stress.piola_kirchoff_2(deformation)
                       for stress, deformation in zip(stresses, deformations)]

        d['fitting_data'] = {'cauchy_stresses': stresses,
                             'eq_stress': eq_stress,
                             'strains': strains,
                             'pk_stresses': pk_stresses,
                             'deformations': deformations
                             }

        logger.info("Analyzing stress/strain data")
        # TODO: @montoyjh: what if it's a cubic system? don't need 6. -computron
        # TODO: Can add population method but want to think about how it should
        #           be done. -montoyjh
        order = self.get('order', 2)
        if order > 2:
            method = 'finite_difference'
        else:
            method = self.get('fitting_method', 'finite_difference')

        if method == 'finite_difference':
            result = ElasticTensorExpansion.from_diff_fit(
                    strains, pk_stresses, eq_stress=eq_stress, order=order)
            if order == 2:
                result = ElasticTensor(result[0])
        elif method == 'pseudoinverse':
            result = ElasticTensor.from_pseudoinverse(strains, pk_stresses)
        elif method == 'independent':
            result = ElasticTensor.from_independent_strains(strains, pk_stresses, eq_stress=eq_stress)
        else:
            raise ValueError("Unsupported method, method must be finite_difference, "
                             "pseudoinverse, or independent")

        ieee = result.convert_to_ieee(ref_struct)
        d.update({
            "elastic_tensor": {
                "raw": result.voigt,
                "ieee_format": ieee.voigt
            }
        })
        if order == 2:
            d.update({"derived_properties": ieee.get_structure_property_dict(ref_struct)})
        else:
            soec = ElasticTensor(ieee[0])
            d.update({"derived_properties": soec.get_structure_property_dict(ref_struct)})

        d["formula_pretty"] = ref_struct.composition.reduced_formula
        d["fitting_method"] = method
        d["order"] = order

        d = jsanitize(d)

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
             "formula_pretty": structure.composition.reduced_formula,
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
                logger.warning("Mode: {} is UNSTABLE. Freq(cm^-1) = {}".format(k, -omega))
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
        poisson (float): poisson ratio. Defaults to 0.25.
        anharmonic_contribution (bool): consider anharmonic contributions to
            Gibbs energy from the Debye model. Defaults to False.
        pressure (float): in GPa, optional.
        metadata (dict): meta data

    """

    required_params = ["tag", "db_file"]
    optional_params = ["qha_type", "t_min", "t_step", "t_max", "mesh", "eos",
                       "pressure", "poisson", "anharmonic_contribution", "metadata"]

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
        anharmonic_contribution = self.get("anharmonic_contribution", False)
        gibbs_dict["metadata"] = self.get("metadata", {})

        db_file = env_chk(self.get("db_file"), fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        # get the optimized structure
        d = mmdb.collection.find_one({"task_label": "{} structure optimization".format(tag)},
                                     {"calcs_reversed": 1})
        structure = Structure.from_dict(d["calcs_reversed"][-1]["output"]['structure'])
        gibbs_dict["structure"] = structure.as_dict()
        gibbs_dict["formula_pretty"] = structure.composition.reduced_formula

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
                                                eos, pressure=pressure, poisson=poisson,
                                                anharmonic_contribution=anharmonic_contribution)
                gibbs_dict.update(qhda.get_summary_dict())
                gibbs_dict["anharmonic_contribution"] = anharmonic_contribution
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

            logger.warning("Quasi-harmonic analysis failed!")
            gibbs_dict["success"] = False
            gibbs_dict["traceback"] = traceback.format_exc()
            gibbs_dict['metadata'].update({"task_label_tag": tag})
            gibbs_dict["created_at"] = datetime.utcnow()

        gibbs_dict = jsanitize(gibbs_dict)

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
        summary_dict["formula_pretty"] = structure.composition.reduced_formula

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
        summary_dict["formula_pretty"] = structure.composition.reduced_formula

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


@explicit_serialize
class MagneticOrderingsToDb(FiretaskBase):
    """
    Used to aggregate tasks docs from magnetic ordering workflow.
    For large-scale/high-throughput use, would recommend a specific
    builder, this is intended for easy, automated use for calculating
    magnetic orderings directly from the get_wf_magnetic_orderings
    workflow. It's unlikely you will want to call this directly.
    
    Required parameters:
        db_file (str): path to the db file that holds your tasks
            collection and that you want to hold the magnetic_orderings
            collection
        wf_uuid (str): auto-generated from get_wf_magnetic_orderings,
            used to make it easier to retrieve task docs
        parent_structure (Structure): Structure of parent crystal (not
            magnetically ordered)
        perform_bader (bool): Perform Bader charge analysis.
        scan (bool): Do static calcs with SCAN functional.

    Optional parameters:
        origins (list): str indicating transformations that generated
            orderings.
        input_index (int): index of input structure to enumerator.
        to_db (bool): if True, the data will be inserted into
            dedicated collection in database, otherwise, will be dumped
            to a .json file.
        additional_fields (dict): fields added to the document such as
            user-defined tags or name, ids, etc

    """

    required_params = ["db_file", "wf_uuid", "parent_structure",
                       "perform_bader", "scan"]
    optional_params = ["origins", "input_index", "to_db", "additional_fields"]

    def run_task(self, fw_spec):
        additional_fields = self.get("additional_fields", {})

        uuid = self["wf_uuid"]
        db_file = env_chk(self.get("db_file"), fw_spec)
        to_db = self.get("to_db", True)

        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)

        formula = self["parent_structure"].formula
        formula_pretty = self["parent_structure"].composition.reduced_formula

        # get ground state energy
        task_label_regex = 'static' if not self['scan'] else 'optimize'
        docs = list(mmdb.collection.find({"wf_meta.wf_uuid": uuid,
                                          "task_label": {"$regex": task_label_regex}},
                                         ["task_id", "output.energy_per_atom"]))

        energies = [d["output"]["energy_per_atom"] for d in docs]
        ground_state_energy = min(energies)
        idx = energies.index(ground_state_energy)
        ground_state_task_id = docs[idx]["task_id"]
        if energies.count(ground_state_energy) > 1:
            logger.warning("Multiple identical energies exist, "
                        "duplicate calculations for {}?".format(formula))

        # get results for different orderings
        docs = list(mmdb.collection.find({
            "task_label": {"$regex": task_label_regex},
            "wf_meta.wf_uuid": uuid
        }))

        summaries = []

        for d in docs:

            # Check if optimizations were done
            if additional_fields.get("relax", True):
                optimize_task_label = d["task_label"].replace("static", "optimize")
                optimize_task = dict(mmdb.collection.find_one({
                            "wf_meta.wf_uuid": uuid,
                            "task_label": optimize_task_label
                        }))
                # used to determine if ordering changed during relaxation
                original_task = optimize_task
                # stored for checking suitable convergence is reached
                energy_diff_relax_static = optimize_task["output"]["energy_per_atom"] \
                                                       - d["output"]["energy_per_atom"]
            else:
                original_task = d
                energy_diff_relax_static = None

            input_structure = Structure.from_dict(optimize_task['input']['structure'])
            input_magmoms = optimize_task['input']['incar']['MAGMOM']
            input_structure.add_site_property('magmom', input_magmoms)

            final_structure = Structure.from_dict(d["output"]["structure"])

            # picking a fairly large threshold so that default 0.6 µB magmoms don't
            # cause problems with analysis, this is obviously not approriate for
            # some magnetic structures with small magnetic moments (e.g. CuO)
            input_analyzer = CollinearMagneticStructureAnalyzer(input_structure, threshold=0.61)
            final_analyzer = CollinearMagneticStructureAnalyzer(final_structure, threshold=0.61)

            if d["task_id"] == ground_state_task_id:
                stable = True
                decomposes_to = None
            else:
                stable = False
                decomposes_to = ground_state_task_id
            energy_above_ground_state_per_atom = d["output"]["energy_per_atom"] \
                                                 - ground_state_energy

            # tells us the order in which structure was guessed
            # 1 is FM, then AFM..., -1 means it was entered manually
            # useful to give us statistics about how many orderings
            # we actually need to calculate
            task_label = d["task_label"].split(' ')
            ordering_index = task_label.index('ordering')
            ordering_index = int(task_label[ordering_index + 1])
            if self.get("origins", None):
                ordering_origin = self["origins"][ordering_index]
            else:
                ordering_origin = None

            final_magmoms = final_structure.site_properties["magmom"]
            magmoms = {"vasp": final_magmoms}
            if self["perform_bader"]:
                # if bader has already been run during task ingestion,
                # use existing analysis
                if "bader" in d:
                    magmoms["bader"] = d["bader"]["magmom"]
                # else try to run it
                else:
                    try:
                        dir_name = d["dir_name"]
                        # strip hostname if present, implicitly assumes
                        # ToDB task has access to appropriate dir
                        if ":" in dir_name:
                            dir_name = dir_name.split(":")[1]
                        magmoms["bader"] = bader_analysis_from_path(dir_name)["magmom"]
                        # prefer bader magmoms if we have them
                        final_magmoms = magmoms["bader"]
                    except Exception as e:
                        magmoms["bader"] = "Bader analysis failed: {}".format(e)

            input_order_check = [0 if abs(m) < 0.61 else m for m in input_magmoms]
            final_order_check = [0 if abs(m) < 0.61 else m for m in final_magmoms]
            ordering_changed = not np.array_equal(np.sign(input_order_check),
                                                  np.sign(final_order_check))

            symmetry_changed = (final_structure.get_space_group_info()[0]
                                != input_structure.get_space_group_info()[0])

            total_magnetization = abs(d["calcs_reversed"][0]["output"]["outcar"]["total_magnetization"])
            num_formula_units = sum(d["calcs_reversed"][0]["composition_reduced"].values())/\
                                sum(d["calcs_reversed"][0]["composition_unit_cell"].values())
            total_magnetization_per_formula_unit = total_magnetization/num_formula_units
            total_magnetization_per_unit_volume = total_magnetization/final_structure.volume

            summary = {
                "formula": formula,
                "formula_pretty": formula_pretty,
                "parent_structure": self["parent_structure"].as_dict(),
                "wf_meta": d["wf_meta"],  # book-keeping
                "task_id": d["task_id"],
                "structure": final_structure.as_dict(),
                "magmoms": magmoms,
                "input": {
                    "structure": input_structure.as_dict(),
                    "ordering": input_analyzer.ordering.value,
                    "symmetry": input_structure.get_space_group_info()[0],
                    "index": ordering_index,
                    "origin": ordering_origin,
                    "input_index": self.get("input_index", None)
                },
                "total_magnetization": total_magnetization,
                "total_magnetization_per_formula_unit": total_magnetization_per_formula_unit,
                "total_magnetization_per_unit_volume": total_magnetization_per_unit_volume,
                "ordering": final_analyzer.ordering.value,
                "ordering_changed": ordering_changed,
                "symmetry": final_structure.get_space_group_info()[0],
                "symmetry_changed": symmetry_changed,
                "energy_per_atom": d["output"]["energy_per_atom"],
                "stable": stable,
                "decomposes_to": decomposes_to,
                "energy_above_ground_state_per_atom": energy_above_ground_state_per_atom,
                "energy_diff_relax_static": energy_diff_relax_static,
                "created_at": datetime.utcnow()
            }

            if fw_spec.get("tags", None):
                summary["tags"] = fw_spec["tags"]

            summary["additional_fields"] = additional_fields

            summaries.append(summary)

        mmdb.collection = mmdb.db["magnetic_orderings"]
        mmdb.collection.insert(summaries)

        logger.info("Magnetic orderings calculation complete.")


@explicit_serialize
class MagneticDeformationToDb(FiretaskBase):
    """
    Used to calculate magnetic deformation from
    get_wf_magnetic_deformation workflow. See docstring
    for that workflow for more information.

    Required parameters:
        db_file (str): path to the db file that holds your tasks
            collection and that you want to hold the magnetic_orderings
            collection
        wf_uuid (str): auto-generated from get_wf_magnetic_orderings,
            used to make it easier to retrieve task docs

    Optional parameters:
        to_db (bool): if True, the data will be inserted into
            dedicated collection in database, otherwise, will be dumped
            to a .json file.

    """

    required_params = ["db_file", "wf_uuid"]
    optional_params = ["to_db"]

    def run_task(self, fw_spec):

        uuid = self["wf_uuid"]
        db_file = env_chk(self.get("db_file"), fw_spec)
        to_db = self.get("to_db", True)

        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)

        # get the non-magnetic structure
        d_nm = mmdb.collection.find_one({
            "task_label": "magnetic deformation optimize non-magnetic",
            "wf_meta.wf_uuid": uuid
        })
        nm_structure = Structure.from_dict(d_nm["output"]["structure"])
        nm_run_stats = d_nm["run_stats"]["overall"]

        # get the magnetic structure
        d_m = mmdb.collection.find_one({
            "task_label": "magnetic deformation optimize magnetic",
            "wf_meta.wf_uuid": uuid
        })
        m_structure = Structure.from_dict(d_m["output"]["structure"])
        m_run_stats = d_m["run_stats"]["overall"]

        msa = CollinearMagneticStructureAnalyzer(m_structure)
        success = False if msa.ordering == Ordering.NM else True

        # calculate magnetic deformation
        mag_def = magnetic_deformation(nm_structure, m_structure).deformation

        # get run stats (mostly used for benchmarking)
        # using same approach as VaspDrone
        try:
            run_stats = {'nm': nm_run_stats, 'm': m_run_stats}
            overall_run_stats = {}
            for key in ["Total CPU time used (sec)", "User time (sec)", "System time (sec)",
                        "Elapsed time (sec)"]:
                overall_run_stats[key] = sum([v[key] for v in run_stats.values()])
        except:
            logger.error("Bad run stats for {}.".format(uuid))
            overall_run_stats = "Bad run stats"

        summary = {
            "formula": nm_structure.composition.reduced_formula,
            "success": success,
            "magnetic_deformation": mag_def,
            "non_magnetic_task_id": d_nm["task_id"],
            "non_magnetic_structure": nm_structure.as_dict(),
            "magnetic_task_id": d_m["task_id"],
            "magnetic_structure": m_structure.as_dict(),
            "run_stats": overall_run_stats,
            "created_at": datetime.utcnow()
        }

        if fw_spec.get("tags", None):
            summary["tags"] = fw_spec["tags"]

        # db_file itself is required but the user can choose to pass the results to db or not
        if to_db:
            mmdb.collection = mmdb.db["magnetic_deformation"]
            mmdb.collection.insert_one(summary)
        else:
            with open("magnetic_deformation.json", "w") as f:
                f.write(json.dumps(summary, default=DATETIME_HANDLER))

        logger.info("Magnetic deformation calculation complete.")


@explicit_serialize
class PolarizationToDb(FiretaskBase):
    """
    Recovers the same branch polarization and spontaneous polarization
    for a ferroelectric workflow.
    """

    optional_params = ["db_file"]

    def run_task(self, fw_spec):

        wfid = list(filter(lambda x: 'wfid' in x, fw_spec['tags'])).pop()
        db_file = env_chk(self.get("db_file"), fw_spec)
        vaspdb = VaspCalcDb.from_db_file(db_file, admin=True)

        # ferroelectric workflow groups calculations by generated wfid tag
        polarization_tasks = vaspdb.collection.find({"tags": wfid, "task_label": {"$regex": ".*polarization"}})

        tasks = []
        outcars = []
        structure_dicts = []
        sort_weight = []
        energies_per_atom = []
        energies = []
        zval_dicts = []

        for p in polarization_tasks:
            # Grab data from each polarization task
            energies_per_atom.append(p['calcs_reversed'][0]['output']['energy_per_atom'])
            energies.append(p['calcs_reversed'][0]['output']['energy'])
            tasks.append(p['task_label'])
            outcars.append(p['calcs_reversed'][0]['output']['outcar'])
            structure_dicts.append(p['calcs_reversed'][0]['input']['structure'])
            zval_dicts.append(p['calcs_reversed'][0]['output']['outcar']['zval_dict'])

            # Add weight for sorting
            # Want polarization calculations in order of nonpolar to polar for Polarization object

            # This number needs to be bigger than the number of calculations
            max_sort_weight = 1000000

            if 'nonpolar_polarization' in p['task_label']:
                sort_weight.append(0)
            elif "polar_polarization" in p['task_label']:
                sort_weight.append(max_sort_weight)
            elif "interpolation_" in p['task_label']:
                num = 0
                part = re.findall(r'interpolation_[0-9]+_polarization', p['task_label'])
                if part != []:
                    part2 = re.findall(r'[0-9]+', part.pop())
                    if part2 != []:
                        num = part2.pop()
                sort_weight.append(max_sort_weight - int(num))

        # Sort polarization tasks
        # nonpolar -> interpolation_n -> interpolation_n-1 -> ...  -> interpolation_1 -> polar
        data = zip(tasks, structure_dicts, outcars, energies_per_atom, energies, sort_weight)
        data = sorted(data,key=lambda x: x[-1])

        # Get the tasks, structures, etc in sorted order from the zipped data.
        tasks, structure_dicts, outcars, energies_per_atom, energies, sort_weight = zip(*data)

        structures = [Structure.from_dict(structure) for structure in structure_dicts]

        # If LCALCPOL = True then Outcar will parse and store the pseudopotential zvals.
        zval_dict = zval_dicts.pop()

        # Assumes that we want to calculate the ionic contribution to the dipole moment.
        # VASP's ionic contribution is sometimes strange.
        # See pymatgen.analysis.ferroelectricity.polarization.Polarization for details.
        p_elecs = [outcar['p_elec'] for outcar in outcars]
        p_ions = [get_total_ionic_dipole(structure, zval_dict) for structure in structures]

        polarization = Polarization(p_elecs, p_ions, structures)

        p_change = np.ravel(polarization.get_polarization_change()).tolist()
        p_norm = polarization.get_polarization_change_norm()
        polarization_max_spline_jumps = polarization.max_spline_jumps()
        same_branch = polarization.get_same_branch_polarization_data(convert_to_muC_per_cm2=True)
        raw_elecs, raw_ions = polarization.get_pelecs_and_pions()
        quanta = polarization.get_lattice_quanta(convert_to_muC_per_cm2=True)

        energy_trend = EnergyTrend(energies_per_atom)
        energy_max_spline_jumps = energy_trend.max_spline_jump()

        polarization_dict = {}

        def split_abc(var, var_name):
            d = {}
            for i, j in enumerate('abc'):
                d.update({var_name + "_{}".format(j): np.ravel(var[:, i]).tolist()})
            return d

        # Add some sort of id for the structures? Like cid but more general?
        # polarization_dict.update({'cid': cid})

        # General information
        polarization_dict.update({'pretty_formula': structures[0].composition.reduced_formula})
        polarization_dict.update({'wfid': wfid})
        polarization_dict.update({'task_label_order': tasks})

        # Polarization information
        polarization_dict.update({'polarization_change': p_change})
        polarization_dict.update({'polarization_change_norm': p_norm})
        polarization_dict.update({'polarization_max_spline_jumps': polarization_max_spline_jumps})
        polarization_dict.update(split_abc(same_branch, "same_branch_polarization"))
        polarization_dict.update(split_abc(raw_elecs, "raw_electron_polarization"))
        polarization_dict.update(split_abc(raw_ions, "raw_electron_polarization"))
        polarization_dict.update(split_abc(quanta, "polarization_quanta"))
        polarization_dict.update({"zval_dict": zval_dict})

        # Energy information
        polarization_dict.update({'energy_per_atom_max_spline_jumps': energy_max_spline_jumps})
        polarization_dict.update({"energies": energies})
        polarization_dict.update({"energies_per_atom": energies_per_atom})
        polarization_dict.update({'outcars': outcars})
        polarization_dict.update({"structures": structure_dicts})

        # Write all the info to db.
        coll = vaspdb.db["polarization_tasks"]
        coll.insert_one(polarization_dict)


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
