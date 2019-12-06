# coding: utf-8


# This module defines tasks that support running QChem in various ways.


import shutil
import os
import subprocess
import logging
import warnings

from pymatgen.io.qchem.inputs import QCInput
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
__date__ = "11/20/18"


logger = get_logger(__name__)

@explicit_serialize
class RunCritic2(FiretaskBase):
    """


    """

    required_params = ["molecule", "cube_file"]

    def run_task(self, fw_spec):
        molecule = self.get("molecule")
        cube = self.get("cube_file")

        compress_at_end = False

        if cube[-3:] == ".gz":
            compress_at_end = True
            decompress_file(cube)
            cube = cube[:-3]

        input_script = ["molecule "+cube]
        input_script += ["load "+cube]
        input_script += ["auto"]
        input_script += ["CPREPORT CP.json"]
        input_script += ["YT JSON YT.json"]
        input_script += ["end"]
        input_script += [""]
        input_script = "\n".join(input_script)

        with open('input_script.cri', 'w') as f:
            f.write(input_script)
        args = ["critic2", "input_script.cri"]

        rs = subprocess.Popen(args,
                              stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              close_fds=True)

        stdout, stderr = rs.communicate()
        stdout = stdout.decode()

        if stderr:
            stderr = stderr.decode()
            warnings.warn(stderr)
            with open('stdout.cri', 'w') as f:
                f.write(stdout)
            with open('stderr.cri', 'w') as f:
                f.write(stderr)

        if rs.returncode != 0:
            raise RuntimeError("critic2 exited with return code {}.".format(rs.returncode))

        CP = loadfn("CP.json")
        bohr_to_ang = 0.529177249

        species = {}
        for specie in CP["structure"]["species"]:
            if specie["name"][1] == "_":
                species[specie["id"]] = specie["name"][0]
            else:
                species[specie["id"]] = specie["name"]

        atoms = []
        centering_vector = CP["structure"]["molecule_centering_vector"]
        for ii,atom in enumerate(CP["structure"]["nonequivalent_atoms"]):
            specie = species[atom["species"]]
            atoms.append(specie)
            tmp = atom["cartesian_coordinates"]
            coords = []
            for jj,val in enumerate(tmp):
                coords.append((val+centering_vector[jj])*bohr_to_ang)
            if str(molecule[ii].specie) != specie:
                raise RuntimeError("Atom ordering different!")
            if molecule[ii].distance_from_point(coords) > 1*10**-6:
                raise RuntimeError("Atom position "+str(ii)+" inconsistent!")

        assert CP["critical_points"]["number_of_nonequivalent_cps"] == CP["critical_points"]["number_of_cell_cps"]

        bond_dict = {}
        for cp in CP["critical_points"]["nonequivalent_cps"]:
            if cp["rank"] == 3 and cp["signature"] == -1:
                bond_dict[cp["id"]] = {"field":cp["field"]}

        for cp in CP["critical_points"]["cell_cps"]:
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
            if bond_dict[cpid]["field"] > 0.02 and bond_dict[cpid]["distance"] < 2.5:
                bonds.append([int(entry)-1 for entry in bond_dict[cpid]["atom_ids"]])

        YT = loadfn("YT.json")
        charges = []
        for site in YT["integration"]["attractors"]:
            charges.append(site["atomic_number"]-site["integrals"][0])

        processed_dict = {}
        processed_dict["bonds"] = bonds
        processed_dict["charges"] = charges
        dumpfn(processed_dict,"processed_critic2.json")

        if compress_at_end:
            compress_file(cube)


