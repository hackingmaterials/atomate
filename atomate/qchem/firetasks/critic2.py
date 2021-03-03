# coding: utf-8


# This module defines a Firetask that runs Critic2 to analyze a Q-Chem electron density.


import shutil
import os
import subprocess
import logging
import warnings

from pymatgen.io.qchem.inputs import QCInput
from pymatgen.command_line.critic2_caller import Critic2Caller
from monty.serialization import loadfn, dumpfn

from custodian import Custodian
from custodian.qchem.handlers import QChemErrorHandler
from custodian.qchem.jobs import QCJob

from monty.tempfile import ScratchDir
from monty.shutil import compress_file, decompress_file
from fireworks import explicit_serialize, FiretaskBase

from atomate.utils.utils import env_chk, get_logger
import numpy as np

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "11/20/19"


logger = get_logger(__name__)

@explicit_serialize
class RunCritic2(FiretaskBase):
    """
    Run the Critic2 package on an electron density cube file produced by a Q-Chem single point calculation
    to generate CP and YT files for electron density critical points analysis.

    Required params:
        molecule (Molecule): Molecule object of the molecule whose electron density is being analyzed
                             Note that if prev_calc_molecule is set in the firework spec it will override
                             the molecule required param.
        cube_file (str): Name of the cube file being analyzed

    """

    required_params = ["molecule", "cube_file"]

    def run_task(self, fw_spec):
        if fw_spec.get("prev_calc_molecule"):
            molecule = fw_spec.get("prev_calc_molecule")
        else:
            molecule = self.get("molecule")
        if molecule == None:
            raise ValueError("No molecule passed and no prev_calc_molecule found in spec! Exiting...")

        compress_at_end = False

        cube = self.get("cube_file")

        if cube[-3:] == ".gz":
            compress_at_end = True
            decompress_file(cube)
            cube = cube[:-3]

        input_script = ["molecule " + cube]
        input_script += ["load " + cube]
        input_script += ["auto"]
        input_script += ["CPREPORT cpreport.json"]
        input_script += ["YT JSON yt.json"]
        input_script += ["end"]
        input_script += [""]
        input_script = "\n".join(input_script)

        caller = Critic2Caller(input_script)

        if compress_at_end:
            compress_file(cube)


@explicit_serialize
class ProcessCritic2(FiretaskBase):
    """
    Process the CP and YT json outputs from a Critic2 execution

    Required params:
        molecule (Molecule): Molecule object of the molecule whose electron density is being analyzed
                             Note that if prev_calc_molecule is set in the firework spec it will override
                             the molecule required param.
    """

    required_params = ["molecule"]

    def run_task(self, fw_spec):
        if fw_spec.get("prev_calc_molecule"):
            molecule = fw_spec.get("prev_calc_molecule")
        else:
            molecule = self.get("molecule")
        if molecule == None:
            raise ValueError("No molecule passed and no prev_calc_molecule found in spec! Exiting...")

        cp_loaded = loadfn("cpreport.json")
        bohr_to_ang = 0.529177249

        species = {}
        for specie in cp_loaded["structure"]["species"]:
            if specie["name"][1] == "_":
                species[specie["id"]] = specie["name"][0]
            else:
                species[specie["id"]] = specie["name"]

        atoms = []
        centering_vector = cp_loaded["structure"]["molecule_centering_vector"]
        for ii,atom in enumerate(cp_loaded["structure"]["nonequivalent_atoms"]):
            specie = species[atom["species"]]
            atoms.append(specie)
            tmp = atom["cartesian_coordinates"]
            coords = []
            for jj,val in enumerate(tmp):
                coords.append((val+centering_vector[jj])*bohr_to_ang)
            if str(molecule[ii].specie) != specie:
                raise RuntimeError("Atom ordering different!")
            if molecule[ii].distance_from_point(coords) > 1*10**-5:
                raise RuntimeError("Atom position "+str(ii)+" inconsistent!")

        if (
            cp_loaded["critical_points"]["number_of_nonequivalent_cps"] !=
            cp_loaded["critical_points"]["number_of_cell_cps"]
        ):
            raise ValueError("ERROR: number_of_nonequivalent_cps should always equal number_of_cell_cps!")

        bond_dict = {}
        for cp in cp_loaded["critical_points"]["nonequivalent_cps"]:
            if cp["rank"] == 3 and cp["signature"] == -1:
                bond_dict[cp["id"]] = {"field":cp["field"]}

        for cp in cp_loaded["critical_points"]["cell_cps"]:
            if cp["id"] in bond_dict:
                # Check if any bonds include fictitious atoms
                bad_bond = False
                for entry in cp["attractors"]:
                    if int(entry["cell_id"])-1 >= len(atoms):
                        bad_bond = True
                # If so, remove them from the bond_dict
                if bad_bond:
                    bond_dict.pop(cp["id"])
                else:
                    bond_dict[cp["id"]]["atom_ids"] = [entry["cell_id"] for entry in cp["attractors"]]
                    bond_dict[cp["id"]]["atoms"] = [atoms[int(entry["cell_id"])-1] for entry in cp["attractors"]]
                    bond_dict[cp["id"]]["distance"] = cp["attractors"][0]["distance"]*bohr_to_ang+cp["attractors"][1]["distance"]*bohr_to_ang
        dumpfn(bond_dict,"bonding.json")

        bonds = []
        for cpid in bond_dict:
            # identify and throw out fictitious bonds
            # NOTE: this should be re-examined and refined in the future
            if bond_dict[cpid]["atoms"] == ["Li","C"] or bond_dict[cpid]["atoms"] == ["C","Li"]:
                if bond_dict[cpid]["field"] > 0.012 and bond_dict[cpid]["distance"] < 2.5:
                    bonds.append([int(entry)-1 for entry in bond_dict[cpid]["atom_ids"]])
            elif bond_dict[cpid]["field"] > 0.02 and bond_dict[cpid]["distance"] < 2.5:
                bonds.append([int(entry)-1 for entry in bond_dict[cpid]["atom_ids"]])

        yt = loadfn("yt.json")
        charges = []
        for site in yt["integration"]["attractors"]:
            charges.append(site["atomic_number"]-site["integrals"][0])

        processed_dict = {}
        processed_dict["bonds"] = bonds
        processed_dict["charges"] = charges
        dumpfn(processed_dict,"processed_critic2.json")
