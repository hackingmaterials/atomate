# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

# This module defines a task that returns all fragments of a molecule

import copy
from pymatgen.core.structure import Molecule
from pymatgen.analysis.ion_placer import IonPlacer
from atomate.utils.utils import env_chk
from atomate.qchem.powerups import use_fake_qchem
from fireworks import FiretaskBase, FWAction, explicit_serialize

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "11/7/18"


@explicit_serialize
class PlaceIon(FiretaskBase):
    """
    """

    optional_params = [
        "molecule", "mulliken", "ion", "charges", "stop_num", "do_triplets", "linked", "qchem_input_params", "test_positions", "ref_dirs"
    ]

    def run_task(self, fw_spec):
        if fw_spec.get("prev_calc_molecule"):
            self.mol = fw_spec.get("prev_calc_molecule")
        elif self.get("molecule"):
            self.mol = self.get("molecule")
        else:
            raise KeyError(
                "No molecule present, add as an optional param or check fw_spec"
            )

        mulliken = None
        if fw_spec.get("prev_calc_mulliken"):
            mulliken = fw_spec.get("prev_calc_mulliken")
        elif self.get("mulliken"):
            mulliken = self.get("mulliken")

        if mulliken == None:
            for site in self.mol:
                if "charge" not in site.properties:
                    raise KeyError("If mulliken not set, each site in the input molecule must already have the charge property! Exiting...")
        elif self.mol.spin_multiplicity != 1:
            self.mol.add_site_property("charge",mulliken[0][::,0])
            self.mol.add_site_property("spin",mulliken[0][::,1])
        else:
            self.mol.add_site_property("charge",mulliken[0])


        self.charges = self.get("charges", [0])
        self.ion = self.get("ion", "Li")
        self.do_triplets = self.get("do_triplets", True)
        self.linked = self.get("linked", False)
        self.qchem_input_params = self.get("qchem_input_params", {})
        if self.get("test_positions") and self.get("ref_dirs"):
            self.testing = True
            self.ion_positions = self.get("test_positions")
        else:
            self.testing = False
            self.ion_positions = IonPlacer(molecule=self.mol, ion=self.ion, stop_num=self.get("stop_num", 100000)).accepted_points
        print(self.ion_positions)
        self._build_molecules()

        return FWAction(detours=self._build_new_FWs())

    def _build_molecules(self):
        """
        Build molecule objects given ion positions, charges, and valid multiplicites
        """
        self.new_molecules = []
        for point in self.ion_positions:
            mol = copy.deepcopy(self.mol)
            mol.remove_site_property("charge")
            mol.remove_site_property("reactivity")
            if mol.spin_multiplicity != 1:
                mol.remove_site_property("spin")
            mol.append(self.ion, point)
            for charge in self.charges:
                new_mol = copy.deepcopy(mol)
                new_mol.set_charge_and_spin(charge=charge)
                self.new_molecules.append(new_mol)
        if self.do_triplets:
            for mol in self.new_molecules:
                if mol.spin_multiplicity == 1:
                    new_mol = copy.deepcopy(mol)
                    new_mol.set_charge_and_spin(charge=mol.charge, spin_multiplicity=3)
                    self.new_molecules.append(new_mol)

    def _build_new_FWs(self):
        """
        Build the list of new fireworks
        """
        from atomate.qchem.fireworks.core import FrequencyFlatteningOptimizeFW
        new_FWs = []
        for ii, molecule in enumerate(self.new_molecules):
            new_FWs.append(
                FrequencyFlatteningOptimizeFW(
                    molecule=molecule,
                    name="ion_pos_" + str(ii),
                    qchem_cmd=">>qchem_cmd<<",
                    max_cores=">>max_cores<<",
                    qchem_input_params=self.qchem_input_params,
                    linked=self.linked,
                    db_file=">>db_file<<"))
        # Extremely jank, but its very hard to test dynamic workflows:
        if self.testing:
            print("testing!")
            from fireworks import Workflow
            tmp_wf = Workflow(new_FWs, name="tmp")
            fake_wf = use_fake_qchem(tmp_wf,self.get("ref_dirs"),"mol.qin.orig.gz")
            new_FWs = fake_wf.fws
        return new_FWs
