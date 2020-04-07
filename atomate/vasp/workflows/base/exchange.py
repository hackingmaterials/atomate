# coding: utf-8

import os
import numpy as np

from fireworks import Workflow, Firework
from atomate.vasp.fireworks.core import StaticFW
from atomate.vasp.powerups import (
    add_tags,
    add_modify_incar,
    add_additional_fields_to_taskdocs,
    add_wf_metadata,
    add_common_powerups,
)
from atomate.vasp.workflows.base.core import get_wf
from atomate.vasp.fireworks.exchange import HeisenbergModelFW, VampireCallerFW
from atomate.vasp.firetasks.parse_outputs import MagneticOrderingsToDb
from atomate.vasp.config import VASP_CMD, DB_FILE, ADD_WF_METADATA
from atomate.vasp.workflows.presets.scan import wf_scan_opt
from atomate.utils.utils import get_logger

from uuid import uuid4
from copy import deepcopy

logger = get_logger(__name__)

from pymatgen.io.vasp.sets import MPStaticSet
from pymatgen.analysis.structure_matcher import StructureMatcher, ElementComparator
from pymatgen.analysis.magnetism.analyzer import (
    CollinearMagneticStructureAnalyzer,
    MagneticStructureEnumerator,
)

__author__ = "Nathan C. Frey"
__maintainer__ = "Nathan C. Frey"
__email__ = "ncfrey@lbl.gov"
__status__ = "Development"
__date__ = "July 2019"
__version__ = 0.1


class ExchangeWF:
    def __init__(
        self,
        magnetic_structures,
        energies,
        default_magmoms=None,
        db_file=DB_FILE,
        name="Exchange WF",
    ):
        """Workflow for computing exchange parameters.

        This workflow takes a set of magnetic orderings and their energies
        from MagneticOrderingsWF and fits to a classical Heisenberg
        Hamiltonian to compute exchange parameters. The critical temperature
        can then be calculated with Monte Carlo.

        Optionally, only the lowest energy FM and AFM configurations can be
        used to compute the average exchange interaction, <J>, without any
        static calculations.

        Args:
            magnetic_structures (list): Structure objects with the 'magmom'
                site property.
            energies (list): Energies per atom in eV.
            default_magmoms (dict): (optional, defaults provided) dict of 
                magnetic elements to their initial magnetic moments in µB, 
                generally these are chosen to be high-spin since they can 
                relax to a low-spin configuration during a DFT electronic
                configuration.
            db_file (string): Path to database file.
            name (string): Name of workflow.

        """

        self.uuid = str(uuid4())
        self.wf_meta = {
            "wf_uuid": self.uuid,
            "wf_name": self.__class__.__name__,
            "wf_version": __version__,
        }

        # Make sure gs is index 0
        ordered_structures = [
            s for _, s in sorted(zip(energies, magnetic_structures), reverse=False)
        ]
        ordered_energies = sorted(energies, reverse=False)

        self.structures = ordered_structures
        self.energies = ordered_energies

        # Check for magmoms
        for s in magnetic_structures:
            try:
                magmoms = s.site_properties["magmom"]
            except:
                raise RuntimeError(
                    "All structures must have 'magmom' site \
                    property."
                )

    def get_wf(self, num_orderings_hard_limit=16, c=None):
        """Retrieve Fireworks workflow.

        c is an optional dictionary that can contain:
        * heisenberg_settings: 
            cutoff (float): Starting point for nearest neighbor search.
            tol (float): Tolerance for equivalent NN bonds.
        * mc_settings:
            mc_box_size (float): MC simulation box size in nm.
            equil_timesteps (int): Number of MC equilibration moves.
            mc_timesteps (int): Number of MC moves for averaging.
            avg (bool): Compute only <J>.
        * DB_FILE:
            path to db.json.

        Args:
            num_orderings_hard_limit (int): will make sure total number of
                magnetic orderings does not exceed this number even if there
                are extra orderings of equivalent symmetry
            c Optional[dict]: additional config dict described above

        Returns: 
            wf (Workflow): Heisenberg Model + Vampire Monte Carlo.

        TODO:
            * Add static SCAN option (only optimization is available)

        """

        c = c or {"DB_FILE": DB_FILE}

        if "DB_FILE" not in c:
            c["DB_FILE"] = DB_FILE

        heisenberg_settings = c.get("heisenberg_settings", {})

        fws = []

        heisenberg_model_fw = HeisenbergModelFW(
            wf_uuid=self.uuid,
            parent_structure=self.structures[0],
            db_file=c["DB_FILE"],
            heisenberg_settings=heisenberg_settings,
            parents=None,
            structures=self.structures,
            energies=self.energies,
        )

        # Vampire Monte Carlo
        mc_settings = c.get("mc_settings", {})

        vampire_fw = VampireCallerFW(
            wf_uuid=self.uuid,
            parent_structure=self.structures[0],
            parents=[heisenberg_model_fw],
            db_file=c["DB_FILE"],
            mc_settings=mc_settings,
        )

        fws = [heisenberg_model_fw, vampire_fw]

        wf = Workflow(fws)

        # Add metadata
        wf = add_additional_fields_to_taskdocs(wf, {"wf_meta": self.wf_meta})
        formula = self.structures[0].composition.reduced_formula
        wf.name = "{} - Exchange".format(formula)

        return wf
