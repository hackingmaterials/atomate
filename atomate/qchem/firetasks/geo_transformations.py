# coding: utf-8


# This module defines firetasks for modifying molecular geometries

import copy

import numpy as np
from fireworks import FiretaskBase, FWAction, explicit_serialize
from pymatgen.io.babel import BabelMolAdaptor

__author__ = "Brandon Wood"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Brandon Wood"
__email__ = "b.wood@berkeley.edu"
__status__ = "Alpha"
__date__ = "5/20/18"
__credits__ = "Sam Blau, Shyam Dwaraknath"


@explicit_serialize
class RotateTorsion(FiretaskBase):
    """
    Writes QChem Input files from input sets. A dictionary is passed to WriteInputFromIOSet where
    parameters are given as keys in the dictionary.

    required_params:
        atom_indexes (list): This should be a list of the pymatgen molecule indexes of the four atoms
        in torsion angle to be rotated
        angle (float): This is the desired torsion angle in degrees (the value should be between -180 and 180)

    optional_params:
        molecule (Molecule): Specify a pymatgen molecule to be rotated. A molecule is optional because molecules can
        be inherited from previous fireworks
    """

    required_params = ["atom_indexes", "angle"]
    optional_params = ["molecule"]

    def run_task(self, fw_spec):
        if fw_spec.get("prev_calc_molecule"):
            start_mol = fw_spec.get("prev_calc_molecule")
        # if a molecule is being passed through fw_spec
        elif self.get("molecule"):
            start_mol = self.get("molecule")
        else:
            raise KeyError(
                "No molecule present, add as an optional param or check fw_spec"
            )

        babe_mol = BabelMolAdaptor(start_mol).openbabel_mol
        babe_mol.SetTorsion(
            self["atom_indexes"][0],
            self["atom_indexes"][1],
            self["atom_indexes"][2],
            self["atom_indexes"][3],
            (self["angle"] * np.pi / 180.0),
        )
        rotated_mol = BabelMolAdaptor(babe_mol).pymatgen_mol

        # update the fw_spec with the rotated geometry
        update_spec = {"prev_calc_molecule": rotated_mol}

        return FWAction(update_spec=update_spec)


@explicit_serialize
class PerturbGeometry(FiretaskBase):
    """
    Modify a Molecule geometry.
    """

    required_params = list()
    optional_params = ["molecule", "mode", "scale"]

    def run_task(self, fw_spec):
        if self.get("molecule"):
            mol = self.get("molecule")
        elif fw_spec.get("prev_calc_molecule"):
            mol = fw_spec.get("prev_calc_molecule")
        else:
            raise KeyError("No molecule present; add as an optional param or check fw_spec")

        if self.get("mode") is not None:
            mode = self.get("mode")
        else:
            raise KeyError("No mode present; add as an optional param or check fw_spec")

        mol_copy = copy.deepcopy(mol)

        for ii in range(len(mol)):
            vec = np.array(mode[ii])
            mol_copy.translate_sites(indices=[ii], vector=vec * self.get("scale", 1.0))

        update_spec = {"prev_calc_molecule": mol_copy}

        return FWAction(update_spec=update_spec)