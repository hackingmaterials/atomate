# coding: utf-8


"""
This Drone tries to produce a more sensible task dictionary than the default VaspToDbTaskDrone.
Some of the changes are documented in this thread:
https://groups.google.com/forum/#!topic/pymatgen/pQ-emBpeV5U
"""

import os
import re
import datetime
from fnmatch import fnmatch
from collections import OrderedDict
import json
import glob
import traceback
import warnings

from monty.io import zopen
from monty.json import jsanitize
from monty.os.path import which

import numpy as np

from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.core.operations import SymmOp
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import BSVasprun, Vasprun, Outcar, Locpot
from pymatgen.io.vasp.inputs import Poscar, Potcar, Incar, Kpoints
from pymatgen.io.vasp.outputs import Chgcar
from pymatgen.apps.borg.hive import AbstractDrone
from pymatgen.command_line.bader_caller import bader_analysis_from_path

from atomate.utils.utils import get_uri

from atomate.utils.utils import get_logger
from atomate import __version__ as atomate_version
from atomate.vasp.config import STORE_VOLUMETRIC_DATA, STORE_ADDITIONAL_JSON

__author__ = 'Kiran Mathew, Shyue Ping Ong, Shyam Dwaraknath, Anubhav Jain'
__email__ = 'kmathew@lbl.gov'
__date__ = 'Mar 27, 2016'
__version__ = "0.1.0"

logger = get_logger(__name__)

bader_exe_exists = which("bader") or which("bader.exe")


class VaspDrone(AbstractDrone):
    """
    pymatgen-db VaspToDbTaskDrone with updated schema and documents processing methods.
    Please refer to matgendb.creator.VaspToDbTaskDrone documentation.
    """

    __version__ = atomate_version  # note: the version is inserted into the task doc

    # Schema def of important keys and sub-keys; used in validation
    schema = {
        "root": {
            "schema", "dir_name", "chemsys", "composition_reduced",
            "formula_pretty", "formula_reduced_abc", "elements",
            "nelements", "formula_anonymous", "calcs_reversed", "completed_at",
            "nsites", "composition_unit_cell", "input", "output", "state",
            "analysis", "run_stats"
        },
        "input": {'is_lasph', 'is_hubbard', 'xc_override', 'potcar_spec',
                  'hubbards', 'structure', 'pseudo_potential'},
        "output": {'structure', 'spacegroup', 'density', 'energy',
                   'energy_per_atom', 'is_gap_direct', 'bandgap', 'vbm',
                   'cbm', 'is_metal', 'forces', 'stress'},
        "calcs_reversed": {
            'dir_name', 'run_type', 'elements', 'nelements',
            'formula_pretty', 'formula_reduced_abc', 'composition_reduced',
            'vasp_version', 'formula_anonymous', 'nsites',
            'composition_unit_cell', 'completed_at', 'task', 'input', 'output',
            'has_vasp_completed'
        },
        "analysis": {'delta_volume_as_percent', 'delta_volume', 'max_force',
                     'errors', 'warnings'}

    }

    def __init__(self, runs=None, parse_dos="auto", bandstructure_mode="auto",
                 parse_locpot=True, additional_fields=None, use_full_uri=True,
                 parse_bader=bader_exe_exists, parse_chgcar=False, parse_aeccar=False,
                 parse_potcar_file=True,
                 store_volumetric_data=STORE_VOLUMETRIC_DATA,
                 store_additional_json=STORE_ADDITIONAL_JSON):
        """
        Initialize a Vasp drone to parse vasp outputs
        Args:
            runs (list): Naming scheme for multiple calcuations in on folder e.g. ["relax1","relax2"].
             Can be subfolder or extension
            parse_dos (str or bool): Whether to parse the DOS. Can be "auto", True or False.
            "auto" will only parse DOS if NSW = 0, so there are no ionic steps
            bandstructure_mode (str or bool): How to parse the bandstructure or not. Can be "auto","line", True or False.
             "auto" will parse the bandstructure with projections for NSCF calcs and decide automatically
              if it's line mode or uniform. Saves the bandstructure in the output doc.
             "line" will parse the bandstructure as a line mode calculation with projections.
              Saves the bandstructure in the output doc.
             True will parse the bandstructure with projections as a uniform calculation.
              Saves the bandstructure in the output doc.
             False will parse the bandstructure without projections to calculate vbm, cbm, band_gap, is_metal and efermi
              Dose not saves the bandstructure in the output doc.
            parse_locpot (bool): Parses the LOCPOT file and saves the 3 axis averages
            additional_fields (dict): dictionary of additional fields to add to output document
            use_full_uri (bool): converts the directory path to the full URI path
            parse_bader (bool): Run and parse Bader charge data. Defaults to True if Bader is present
            parse_chgcar (bool): Run and parse CHGCAR file
            parse_aeccar (bool): Run and parse AECCAR0 and AECCAR2 files
            store_volumetric_data (list): List of files to store, choose from ('CHGCAR', 'LOCPOT',
            'AECCAR0', 'AECCAR1', 'AECCAR2', 'ELFCAR'), case insensitive
            store_additional_json (bool): If True, parse any .json files present and store as
            sub-doc including the FW.json if present
        """
        self.parse_dos = parse_dos
        self.additional_fields = additional_fields or {}
        self.use_full_uri = use_full_uri
        self.runs = runs or ["precondition"] + ["relax" + str(i + 1)
                                                for i in range(9)]  # can't auto-detect: path unknown
        self.bandstructure_mode = bandstructure_mode
        self.parse_locpot = parse_locpot
        self.parse_bader = parse_bader
        self.store_volumetric_data = [f.lower() for f in store_volumetric_data]
        self.store_additional_json = store_additional_json
        self.parse_potcar_file = parse_potcar_file

        if parse_chgcar or parse_aeccar:
            warnings.warn("These options have been deprecated in favor of the 'store_volumetric_data' "
                          "keyword argument, which is more general. Functionality is equivalent.",
                          DeprecationWarning)
            if parse_chgcar and "chgcar" not in self.store_volumetric_data:
                self.store_volumetric_data.append("chgcar")
            if parse_aeccar and "aeccar0" not in self.store_volumetric_data:
                self.store_volumetric_data.append("aeccar0")
            if parse_aeccar and "aeccar2" not in self.store_volumetric_data:
                self.store_volumetric_data.append("aeccar2")

    def assimilate(self, path):
        """
        Adapted from matgendb.creator
        Parses vasp runs(vasprun.xml file) and insert the result into the db.
        Get the entire task doc from the vasprum.xml and the OUTCAR files in the path.
        Also adds some post-processed info.

        Args:
            path (str): Path to the directory containing vasprun.xml and OUTCAR files

        Returns:
            (dict): a task dictionary
        """
        logger.info("Getting task doc for base dir :{}".format(path))
        vasprun_files = self.filter_files(path, file_pattern="vasprun.xml")
        outcar_files = self.filter_files(path, file_pattern="OUTCAR")
        if len(vasprun_files) > 0 and len(outcar_files) > 0:
            d = self.generate_doc(path, vasprun_files, outcar_files)
            self.post_process(path, d)
        else:
            raise ValueError("No VASP files found!")
        self.validate_doc(d)
        return d

    def filter_files(self, path, file_pattern="vasprun.xml"):
        """
        Find the files that match the pattern in the given path and
        return them in an ordered dictionary. The searched for files are
        filtered by the run types defined in self.runs. e.g. ["relax1", "relax2", ...].
        Only 2 schemes of the file filtering is enabled: searching for run types
        in the list of files and in the filenames. Modify this method if more
        sophisticated filtering scheme is needed.

        Args:
            path (string): path to the folder
            file_pattern (string): files to be searched for

        Returns:
            OrderedDict of the names of the files to be processed further.
            The key is set from list of run types: self.runs
        """
        processed_files = OrderedDict()
        files = os.listdir(path)
        for r in self.runs:
            # try subfolder schema
            if r in files:
                for f in os.listdir(os.path.join(path, r)):
                    if fnmatch(f, "{}*".format(file_pattern)):
                        processed_files[r] = os.path.join(r, f)
            # try extension schema
            else:
                for f in files:
                    if fnmatch(f, "{}.{}*".format(file_pattern, r)):
                        processed_files[r] = f
        if len(processed_files) == 0:
            # get any matching file from the folder
            for f in files:
                if fnmatch(f, "{}*".format(file_pattern)):
                    processed_files['standard'] = f
        return processed_files

    def generate_doc(self, dir_name, vasprun_files, outcar_files):
        """
        Adapted from matgendb.creator.generate_doc
        """
        try:
            # basic properties, incl. calcs_reversed and run_stats
            fullpath = os.path.abspath(dir_name)
            d = jsanitize(self.additional_fields, strict=True)
            d["schema"] = {"code": "atomate", "version": VaspDrone.__version__}
            d["dir_name"] = fullpath
            d["calcs_reversed"] = [self.process_vasprun(dir_name, taskname, filename)
                                   for taskname, filename in vasprun_files.items()]
            outcar_data = [Outcar(os.path.join(dir_name, filename)).as_dict()
                           for taskname, filename in outcar_files.items()]
            run_stats = {}
            for i, d_calc in enumerate(d["calcs_reversed"]):
                run_stats[d_calc["task"]["name"]] = outcar_data[i].pop("run_stats")
                if d_calc.get("output"):
                    d_calc["output"].update({"outcar": outcar_data[i]})
                else:
                    d_calc["output"] = {"outcar": outcar_data[i]}
            try:
                overall_run_stats = {}
                for key in ["Total CPU time used (sec)", "User time (sec)", "System time (sec)",
                            "Elapsed time (sec)"]:
                    overall_run_stats[key] = sum([v[key] for v in run_stats.values()])
                run_stats["overall"] = overall_run_stats
            except:
                logger.error("Bad run stats for {}.".format(fullpath))
            d["run_stats"] = run_stats

            # reverse the calculations data order so newest calc is first
            d["calcs_reversed"].reverse()

            # set root formula/composition keys based on initial and final calcs
            d_calc_init = d["calcs_reversed"][-1]
            d_calc_final = d["calcs_reversed"][0]
            d["chemsys"] = "-".join(sorted(d_calc_final["elements"]))
            comp = Composition(d_calc_final["composition_unit_cell"])
            d["formula_anonymous"] = comp.anonymized_formula
            d["formula_reduced_abc"] = comp.reduced_composition.alphabetical_formula
            for root_key in ["completed_at", "nsites", "composition_unit_cell",
                             "composition_reduced", "formula_pretty", "elements", "nelements"]:
                d[root_key] = d_calc_final[root_key]

            # store the input key based on initial calc
            # store any overrides to the exchange correlation functional
            xc = d_calc_init["input"]["incar"].get("GGA")
            if xc:
                xc = xc.upper()
            p = d_calc_init["input"]["potcar_type"][0].split("_")
            pot_type = p[0]
            functional = "lda" if len(pot_type) == 1 else "_".join(p[1:])
            d["input"] = {"structure": d_calc_init["input"]["structure"],
                          "is_hubbard": d_calc_init.pop("is_hubbard"),
                          "hubbards": d_calc_init.pop("hubbards"),
                          "is_lasph": d_calc_init["input"]["incar"].get("LASPH", False),
                          "potcar_spec": d_calc_init["input"].get("potcar_spec"),
                          "xc_override": xc,
                          "pseudo_potential": {"functional": functional.lower(),
                                               "pot_type": pot_type.lower(),
                                               "labels": d_calc_init["input"]["potcar"]},
                          "parameters": d_calc_init["input"]["parameters"],
                          "incar": d_calc_init["input"]["incar"]
                          }

            # store the output key based on final calc
            d["output"] = {
                "structure": d_calc_final["output"]["structure"],
                "density": d_calc_final.pop("density"),
                "energy": d_calc_final["output"]["energy"],
                "energy_per_atom": d_calc_final["output"]["energy_per_atom"],
                "forces": d_calc_final["output"]["ionic_steps"][-1].get("forces"),
                "stress": d_calc_final["output"]["ionic_steps"][-1].get("stress")}

            # patch calculated magnetic moments into final structure
            if len(d_calc_final["output"]["outcar"]["magnetization"]) != 0:
                magmoms = [m["tot"] for m in d_calc_final["output"]["outcar"]["magnetization"]]
                s = Structure.from_dict(d["output"]["structure"])
                s.add_site_property('magmom', magmoms)
                d["output"]["structure"] = s.as_dict()

            calc = d["calcs_reversed"][0]

            # copy band gap and properties into output
            d["output"].update({"bandgap": calc["output"]["bandgap"],
                                "cbm": calc["output"]["cbm"],
                                "vbm": calc["output"]["vbm"],
                                "is_gap_direct": calc["output"]["is_gap_direct"]})
            try:
                d["output"].update({"is_metal": calc["output"]["is_metal"]})
                if not calc["output"]["is_gap_direct"]:
                    d["output"]["direct_gap"] = calc["output"]["direct_gap"]
                if "transition" in calc["output"]:
                    d["output"]["transition"] = calc["output"]["transition"]

            except Exception:
                if self.bandstructure_mode is True:
                    logger.error(traceback.format_exc())
                    logger.error("Error in " + os.path.abspath(dir_name) + ".\n" + traceback.format_exc())
                    raise

            # Store symmetry information
            sg = SpacegroupAnalyzer(Structure.from_dict(d_calc_final["output"]["structure"]), 0.1)
            if not sg.get_symmetry_dataset():
                sg = SpacegroupAnalyzer(Structure.from_dict(d_calc_final["output"]["structure"]),
                                        1e-3, 1)
            d["output"]["spacegroup"] = {
                "source": "spglib",
                "symbol": sg.get_space_group_symbol(),
                "number": sg.get_space_group_number(),
                "point_group": sg.get_point_group_symbol(),
                "crystal_system": sg.get_crystal_system(),
                "hall": sg.get_hall()}

            # store dieelctric and piezo information
            if d["input"]["parameters"].get("LEPSILON"):
                for k in ['epsilon_static', 'epsilon_static_wolfe', 'epsilon_ionic']:
                    d["output"][k] = d_calc_final["output"][k]
                if SymmOp.inversion() not in sg.get_symmetry_operations():
                    for k in ["piezo_ionic_tensor", "piezo_tensor"]:
                        d["output"][k] = d_calc_final["output"]["outcar"][k]

            d["state"] = "successful" if d_calc["has_vasp_completed"] else "unsuccessful"

            self.set_analysis(d)

            d["last_updated"] = datetime.datetime.utcnow()
            return d

        except Exception:
            logger.error(traceback.format_exc())
            logger.error("Error in " + os.path.abspath(dir_name) + ".\n" + traceback.format_exc())
            raise

    def process_vasprun(self, dir_name, taskname, filename):
        """
        Adapted from matgendb.creator

        Process a vasprun.xml file.
        """
        vasprun_file = os.path.join(dir_name, filename)

        vrun = Vasprun(vasprun_file, parse_potcar_file=self.parse_potcar_file)

        d = vrun.as_dict()

        # rename formula keys
        for k, v in {"formula_pretty": "pretty_formula",
                     "composition_reduced": "reduced_cell_formula",
                     "composition_unit_cell": "unit_cell_formula"}.items():
            d[k] = d.pop(v)

        for k in ["eigenvalues", "projected_eigenvalues"]:  # large storage space breaks some docs
            if k in d["output"]:
                del d["output"][k]

        comp = Composition(d["composition_unit_cell"])
        d["formula_anonymous"] = comp.anonymized_formula
        d["formula_reduced_abc"] = comp.reduced_composition.alphabetical_formula
        d["dir_name"] = os.path.abspath(dir_name)
        d["completed_at"] = str(datetime.datetime.fromtimestamp(os.path.getmtime(vasprun_file)))
        d["density"] = vrun.final_structure.density

        # replace 'crystal' with 'structure'
        d["input"]["structure"] = d["input"].pop("crystal")
        d["output"]["structure"] = d["output"].pop("crystal")
        for k, v in {"energy": "final_energy", "energy_per_atom": "final_energy_per_atom"}.items():
            d["output"][k] = d["output"].pop(v)

        # Process bandstructure and DOS
        if self.bandstructure_mode != False:
            bs = self.process_bandstructure(vrun)
            if bs:
                d["bandstructure"] = bs

        if self.parse_dos != False:
            dos = self.process_dos(vrun)
            if dos:
                d["dos"] = dos

        # Parse electronic information if possible.
        # For certain optimizers this is broken and we don't get an efermi resulting in the bandstructure
        try:
            bs = vrun.get_band_structure()
            bs_gap = bs.get_band_gap()
            d["output"]["vbm"] = bs.get_vbm()["energy"]
            d["output"]["cbm"] = bs.get_cbm()["energy"]
            d["output"]["bandgap"] = bs_gap["energy"]
            d["output"]["is_gap_direct"] = bs_gap["direct"]
            d["output"]["is_metal"] = bs.is_metal()
            if not bs_gap["direct"]:
                d["output"]["direct_gap"] = bs.get_direct_band_gap()
            if isinstance(bs, BandStructureSymmLine):
                d["output"]["transition"] = bs_gap["transition"]

        except Exception:
            logger.warning("Error in parsing bandstructure")
            if vrun.incar["IBRION"] == 1:
                logger.warning("Vasp doesn't properly output efermi for IBRION == 1")
            if self.bandstructure_mode is True:
                logger.error(traceback.format_exc())
                logger.error("Error in " + os.path.abspath(dir_name) + ".\n" + traceback.format_exc())
                raise

        # store run name and location ,e.g. relax1, relax2, etc.
        d["task"] = {"type": taskname, "name": taskname}

        # include output file names
        d["output_file_paths"] = self.process_raw_data(dir_name, taskname=taskname)

        # parse axially averaged locpot
        if "locpot" in d["output_file_paths"] and self.parse_locpot:
            locpot = Locpot.from_file(os.path.join(dir_name, d["output_file_paths"]["locpot"]))
            d["output"]["locpot"] = {i: locpot.get_average_along_axis(i) for i in range(3)}

        if self.store_volumetric_data:
            for file in self.store_volumetric_data:
                if file in d["output_file_paths"]:
                    try:
                        # assume volumetric data is all in CHGCAR format
                        data = Chgcar.from_file(os.path.join(dir_name, d["output_file_paths"][file]))
                        d[file] = data
                    except:
                        raise ValueError("Failed to parse {} at {}.".format(file,
                                                                            d["output_file_paths"][file]))

        # parse force constants
        if hasattr(vrun, "force_constants"):
            d["output"]["force_constants"] = vrun.force_constants.tolist()
            d["output"]["normalmode_eigenvals"] = vrun.normalmode_eigenvals.tolist()
            d["output"]["normalmode_eigenvecs"] = vrun.normalmode_eigenvecs.tolist()

        # perform Bader analysis using Henkelman bader
        if self.parse_bader and "chgcar" in d["output_file_paths"]:
            suffix = '' if taskname == 'standard' else ".{}".format(taskname)
            bader = bader_analysis_from_path(dir_name, suffix=suffix)
            d["bader"] = bader

        return d

    def process_bandstructure(self, vrun):

        vasprun_file = vrun.filename
        # Band structure parsing logic
        if str(self.bandstructure_mode).lower() == "auto":
            # if NSCF calculation
            if vrun.incar.get("ICHARG", 0) > 10:
                bs_vrun = BSVasprun(vasprun_file, parse_projected_eigen=True)
                try:
                    # Try parsing line mode
                    bs = bs_vrun.get_band_structure(line_mode=True)
                except:
                    # Just treat as a regular calculation
                    bs = bs_vrun.get_band_structure()
            # else just regular calculation
            else:
                bs_vrun = BSVasprun(vasprun_file, parse_projected_eigen=False)
                bs = bs_vrun.get_band_structure()

            # only save the bandstructure if not moving ions
            if vrun.incar.get("NSW", 0) <= 1:
                return bs.as_dict()

        # legacy line/True behavior for bandstructure_mode
        elif self.bandstructure_mode:
            bs_vrun = BSVasprun(vasprun_file, parse_projected_eigen=True)
            bs = bs_vrun.get_band_structure(line_mode=(str(self.bandstructure_mode).lower() == "line"))
            return bs.as_dict()

        return None

    def process_dos(self, vrun):
        # parse dos if forced to or auto mode set and  0 ionic steps were performed -> static calculation and not DFPT
        if self.parse_dos == True or (str(self.parse_dos).lower() == "auto" and vrun.incar.get("NSW", 0) < 1):
            try:
                return vrun.complete_dos.as_dict()
            except:
                raise ValueError("No valid dos data exist")

    def process_raw_data(self, dir_name, taskname="standard"):
        """
        It is useful to store what raw data has been calculated
        and exists for easier querying of the taskdoc.

        :param dir_name: directory to search
        :param taskname: taskname, e.g. "relax1"
        :return: dict of files present
        """
        d = {}
        possible_files = ('CHGCAR', 'LOCPOT', 'AECCAR0', 'AECCAR1', 'AECCAR2',
                          'ELFCAR', 'WAVECAR', 'PROCAR', 'OPTIC')
        for f in possible_files:
            files = self.filter_files(dir_name, file_pattern=f)
            if taskname in files:
                d[f.lower()] = files[taskname]
        return d

    @staticmethod
    def set_analysis(d, max_force_threshold=0.5, volume_change_threshold=0.2):
        """
        Adapted from matgendb.creator

        set the 'analysis' key
        """
        initial_vol = d["input"]["structure"]["lattice"]["volume"]
        final_vol = d["output"]["structure"]["lattice"]["volume"]
        delta_vol = final_vol - initial_vol
        percent_delta_vol = 100 * delta_vol / initial_vol
        warning_msgs = []
        error_msgs = []

        # delta volume checks
        if abs(percent_delta_vol) > volume_change_threshold:
            warning_msgs.append("Volume change > {}%".format(volume_change_threshold * 100))

        # max force and valid structure checks
        max_force = None
        calc = d["calcs_reversed"][0]
        if d["state"] == "successful":

            # calculate max forces
            if 'forces' in calc['output']['ionic_steps'][-1]:
                forces = np.array(calc['output']['ionic_steps'][-1]['forces'])
                # account for selective dynamics
                final_structure = Structure.from_dict(calc['output']['structure'])
                sdyn = final_structure.site_properties.get('selective_dynamics')
                if sdyn:
                    forces[np.logical_not(sdyn)] = 0
                max_force = max(np.linalg.norm(forces, axis=1))

            if calc["input"]["parameters"].get("NSW", 0) > 0:

                drift = calc['output']['outcar'].get("drift", [[0, 0, 0]])
                max_drift = max([np.linalg.norm(d) for d in drift])
                ediffg = calc["input"]["parameters"].get("EDIFFG", None)
                if ediffg and float(ediffg) < 0:
                    desired_force_convergence = -float(ediffg)
                else:
                    desired_force_convergence = np.inf
                if max_drift > desired_force_convergence:
                    warning_msgs.append("Drift ({}) > desired force convergence ({}), "
                                        "structure likely not converged to desired accuracy."
                                       .format(drift, desired_force_convergence))

            s = Structure.from_dict(d["output"]["structure"])
            if not s.is_valid():
                error_msgs.append("Bad structure (atoms are too close!)")
                d["state"] = "error"

        d["analysis"] = {"delta_volume": delta_vol,
                         "delta_volume_as_percent": percent_delta_vol,
                         "max_force": max_force,
                         "warnings": warning_msgs,
                         "errors": error_msgs}

    def post_process(self, dir_name, d):
        """
        Post-processing for various files other than the vasprun.xml and OUTCAR.
        Looks for files: transformations.json and custodian.json. Modify this if other
        output files need to be processed.

        Args:
            dir_name:
                The dir_name.
            d:
                Current doc generated.
        """
        logger.info("Post-processing dir:{}".format(dir_name))
        fullpath = os.path.abspath(dir_name)
        # VASP input generated by pymatgen's alchemy has a transformations.json file that tracks
        # the origin of a particular structure. If such a file is found, it is inserted into the
        # task doc as d["transformations"]
        transformations = {}
        filenames = glob.glob(os.path.join(fullpath, "transformations.json*"))
        if len(filenames) >= 1:
            with zopen(filenames[0], "rt") as f:
                transformations = json.load(f)
                try:
                    m = re.match("(\d+)-ICSD", transformations["history"][0]["source"])
                    if m:
                        d["icsd_id"] = int(m.group(1))
                except Exception as ex:
                    logger.warning("Cannot parse ICSD from transformations file.")
                    pass
        else:
            logger.warning("Transformations file does not exist.")

        other_parameters = transformations.get("other_parameters")
        new_tags = None
        if other_parameters:
            # We don't want to leave tags or authors in the
            # transformations file because they'd be copied into
            # every structure generated after this one.
            new_tags = other_parameters.pop("tags", None)
            new_author = other_parameters.pop("author", None)
            if new_author:
                d["author"] = new_author
            if not other_parameters:  # if dict is now empty remove it
                transformations.pop("other_parameters")
        d["transformations"] = transformations

        # Calculations done using custodian has a custodian.json,
        # which tracks the jobs performed and any errors detected and fixed.
        # This is useful for tracking what has actually be done to get a
        # result. If such a file is found, it is inserted into the task doc
        # as d["custodian"]
        filenames = glob.glob(os.path.join(fullpath, "custodian.json*"))
        if len(filenames) >= 1:
            custodian = []
            for fname in filenames:
                with zopen(fname, "rt") as f:
                    custodian.append(json.load(f)[0])
            d["custodian"] = custodian
        # Convert to full uri path.
        if self.use_full_uri:
            d["dir_name"] = get_uri(dir_name)
        if new_tags:
            d["tags"] = new_tags

        # Calculations using custodian generate a *.orig file for the inputs
        # This is useful to know how the calculation originally started
        # if such files are found they are inserted into orig_inputs
        filenames = glob.glob(os.path.join(fullpath, "*.orig*"))

        if len(filenames) >= 1:
            d["orig_inputs"] = {}
            for f in filenames:
                if "INCAR.orig" in f:
                    d["orig_inputs"]["incar"] = Incar.from_file(f).as_dict()
                if "POTCAR.orig" in f:
                    d["orig_inputs"]["potcar"] = Potcar.from_file(f).as_dict()
                if "KPOINTS.orig" in f:
                    d["orig_inputs"]["kpoints"] = Kpoints.from_file(f).as_dict()
                if "POSCAR.orig" in f:
                    d["orig_inputs"]["poscar"] = Poscar.from_file(f).as_dict()

        filenames = glob.glob(os.path.join(fullpath, "*.json*"))
        if self.store_additional_json and filenames:
            for filename in filenames:
                key = os.path.basename(filename).split('.')[0]
                if key != "custodian" and key != "transformations":
                    with zopen(filename, "rt") as f:
                        d[key] = json.load(f)

        logger.info("Post-processed " + fullpath)

    def validate_doc(self, d):
        """
        Sanity check.
        Make sure all the important keys are set
        """
        # TODO: @matk86 - I like the validation but I think no one will notice a failed
        # validation tests which removes the usefulness of this. Any ideas to make people
        # notice if the validation fails? -computron
        for k, v in self.schema.items():
            if k == "calcs_reversed":
                diff = v.difference(set(d.get(k, d)[0].keys()))
            else:
                diff = v.difference(set(d.get(k, d).keys()))
            if diff:
                logger.warning("The keys {0} in {1} not set".format(diff, k))

    def get_valid_paths(self, path):
        """
        There are some restrictions on the valid directory structures:

        1. There can be only one vasp run in each directory. Nested directories
           are fine.
        2. Directories designated "relax1"..."relax9" are considered to be
           parts of a multiple-optimization run.
        3. Directories containing vasp output with ".relax1"...".relax9" are
           also considered as parts of a multiple-optimization run.
        """
        (parent, subdirs, files) = path
        if set(self.runs).intersection(subdirs):
            return [parent]
        if not any([parent.endswith(os.sep + r) for r in self.runs]) and \
                len(glob.glob(os.path.join(parent, "vasprun.xml*"))) > 0:
            return [parent]
        return []

    def as_dict(self):
        init_args = {
            "parse_dos": self.parse_dos,
            "bandstructure_mode": self.bandstructure_mode,
            "additional_fields": self.additional_fields,
            "use_full_uri": self.use_full_uri,
            "runs": self.runs}
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "version": self.__class__.__version__,
                "init_args": init_args
                }

    @classmethod
    def from_dict(cls, d):
        return cls(**d["init_args"])
