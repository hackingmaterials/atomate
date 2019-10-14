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
from atomate.vasp.firetasks.exchange_tasks import ExchangeToDB
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
        vasp_cmd=VASP_CMD,
        db_file=DB_FILE,
        name="Exchange WF",
    ):
        """Workflow for computing exchange parameters.

        This workflow takes a set of magnetic orderings and their energies from MagneticOrderingsWF and runs static calculations of magnetic orderings in the ground state geometry. The static orderings/energies are then fit to a classical Heisenberg Hamiltonian to compute exchange parameters. The critical temperature can then be calculated with Monte Carlo.

        Args:
            magnetic_structures (list): Structure objects with the 'magmom'
        site property.
            energies (list): Energies per atom in eV.
            default_magmoms: (optional, defaults provided) dict of
        magnetic elements to their initial magnetic moments in µB, generally
        these are chosen to be high-spin since they can relax to a low-spin
        configuration during a DFT electronic configuration

        """

        self.uuid = str(uuid4())
        self.wf_meta = {
            "wf_uuid": self.uuid,
            "wf_name": self.__class__.__name__,
            "wf_version": __version__,
        }

        # self.matched_structures = self.magnetic_structures
        self.matched_structures, self.input_index, self.ordered_structure_origins = self.get_commensurate_orderings(
            magnetic_structures, energies
        )

    @staticmethod
    def get_commensurate_orderings(magnetic_structures, energies):
        """Generate supercells for static calculations.

        From the ground state magnetic ordering and first excited state, generate supercells with magnetic orderings.

        If the gs is FM, match supercells to the 1st es.
        If the gs is AFM, FM is included by default, 1st es may not be.
        If the gs is FiM, 1st es may not be captured.

        Args:
            magnetic_structures (list): Structures.
            energies (list): Energies per atom.

        Returns:
            matched_structures (list): Commensurate supercells for static 
        calculations.

        TODO:
            * Only consider orderings with |S_i| = ground state

        """

        # Sort by energies
        ordered_structures = [
            s for _, s in sorted(zip(energies, magnetic_structures), reverse=False)
        ]
        ordered_energies = sorted(energies, reverse=False)

        # Ground state and 1st excited state
        gs_struct = ordered_structures[0]
        es_struct = ordered_structures[1]

        cmsa = CollinearMagneticStructureAnalyzer(
            es_struct, threshold=0.0, make_primitive=False
        )
        es_ordering = cmsa.ordering.value

        cmsa = CollinearMagneticStructureAnalyzer(
            gs_struct, threshold=0.0, make_primitive=False
        )
        gs_ordering = cmsa.ordering.value

        # FM gs will always be commensurate so we match to the 1st es
        if gs_ordering == "FM":
            enum_struct = es_struct
            fm_moments = np.array(gs_struct.site_properties["magmom"])
            fm_moment = np.mean(fm_moments[np.nonzero(fm_moments)])
            gs_magmoms = [fm_moment for m in es_struct.site_properties["magmom"]]
        elif gs_ordering in ["FiM", "AFM"]:
            enum_struct = gs_struct
            gs_magmoms = [abs(m) for m in gs_struct.site_properties["magmom"]]

        mse = MagneticStructureEnumerator(
            enum_struct,
            strategies=("ferromagnetic", "antiferromagnetic"),
            automatic=False,
            transformation_kwargs={
                "min_cell_size": 1,
                "max_cell_size": 2,
                "check_ordered_symmetry": False,
            },
        )

        # Enumerator bookkeeping
        input_index = mse.input_index
        ordered_structure_origins = mse.ordered_structure_origins

        matched_structures = []

        sm = StructureMatcher(
            primitive_cell=False, attempt_supercell=True, comparator=ElementComparator()
        )

        # Get commensurate supercells
        for s in mse.ordered_structures:
            try:
                s2 = sm.get_s2_like_s1(enum_struct, s)
            except:
                s2 = None
            if s2 is not None:
                # Standardize magnetic structure
                cmsa = CollinearMagneticStructureAnalyzer(s2, 
                    threshold=0.0, make_primitive=False)
                s2 = cmsa.structure
                matched_structures.append(s2)

        # Find the gs ordering in the enumerated supercells
        cmsa = CollinearMagneticStructureAnalyzer(
            gs_struct, threshold=0.0, make_primitive=False
        )

        enum_index = [cmsa.matches_ordering(s) for s in matched_structures].index(
            True
        )

        # Enforce all magmom magnitudes to match the gs
        for s in matched_structures:
            ms = s.site_properties["magmom"]
            magmoms = [np.sign(m1) * m2 for m1, m2 in zip(ms, gs_magmoms)]
            s.add_site_property("magmom", magmoms)

        return matched_structures, input_index, ordered_structure_origins

    def get_wf(self, num_orderings_hard_limit=16, c=None):
        """Retrieve Fireworks workflow.

        c is a dictionary that can contain:
        * user_incar_settings
        * heisenberg_settings: 
            cutoff (float): Starting point for nearest neighbor search.
            tol (float): Tolerance for equivalent NN bonds.
        * mc_settings:
            mc_box_size (float): MC simulation box size in nm.
            equil_timesteps (int): Number of MC equilibration moves.
            mc_timesteps (int): Number of MC moves for averaging.

        Args:
            num_orderings_hard_limit (int): will make sure total number of magnetic
        orderings does not exceed this number even if there are extra orderings
        of equivalent symmetry
            c (dict): additional config dict (as used elsewhere in atomate)

        Returns: 
            wf (Workflow): Static magnetic orderings + Heisenberg Model + Vampire Monte Carlo.

        TODO:
            * Add static SCAN option (only optimization is available)
        """

        c = c or {"VASP_CMD": VASP_CMD, "DB_FILE": DB_FILE}

        fws = []
        analysis_parents = []

        matched_structures = self.matched_structures

        matched_structures = matched_structures[:num_orderings_hard_limit]

        # default incar settings
        user_incar_settings = {"ISYM": 0, "LASPH": True}
        user_incar_settings.update(c.get("user_incar_settings", {}))
        c["user_incar_settings"] = user_incar_settings

        static_fws = []
        # Add static vasp calc FWs which go to 'tasks' collection
        for idx, struct in enumerate(matched_structures):

            cmsa = CollinearMagneticStructureAnalyzer(
                struct, threshold=0.0, make_primitive=False
            )

            name = " ordering {} {} -".format(idx, cmsa.ordering.value)

            vis = MPStaticSet(struct, user_incar_settings=user_incar_settings)
            static_fws.append(
                StaticFW(
                    struct,
                    vasp_input_set=vis,
                    vasp_cmd=c["VASP_CMD"],
                    db_file=c["DB_FILE"],
                    name=name + "static",
                )
            )

            analysis_parents.append(static_fws[-1])

        # Magnetic ordering analysis that generates 'magnetic_orderings'
        # collection
        fw_analysis = Firework(
            ExchangeToDB(
                db_file=c["DB_FILE"],
                wf_uuid=self.uuid,
                auto_generated=False,
                name="MagneticOrderingsToDB",
                parent_structure=self.matched_structures[0],
                input_index=self.input_index,
                origins=self.ordered_structure_origins,
                perform_bader=False,
                scan=False,
            ),
            name="Magnetic Exchange Analysis",
            parents=analysis_parents,
            spec={"_allow_fizzled_parents": True},
        )

        # Heisenberg model mapping
        heisenberg_settings = {"cutoff": 3.0, "tol": 0.04}
        heisenberg_settings.update(c.get("heisenberg_settings", {}))
        c["heisenberg_settings"] = heisenberg_settings

        heisenberg_model_fw = HeisenbergModelFW(
            exchange_wf_uuid=self.uuid,
            parent_structure=self.matched_structures[0],
            parents=[fw_analysis],
            db_file=c["DB_FILE"],
            heisenberg_settings=c["heisenberg_settings"],
        )

        # Vampire Monte Carlo
        mc_settings = {
            "mc_box_size": 4.0,
            "equil_timesteps": 2000,
            "mc_timesteps": 4000,
        }
        mc_settings.update(c.get("mc_settings", {}))
        c["mc_settings"] = mc_settings

        vampire_fw = VampireCallerFW(
            exchange_wf_uuid=self.uuid,
            parent_structure=self.matched_structures[0],
            parents=[heisenberg_model_fw],
            db_file=c["DB_FILE"],
            mc_settings=c["mc_settings"],
        )

        # Chain FWs together
        fws = static_fws + [fw_analysis, heisenberg_model_fw, vampire_fw]

        wf = Workflow(fws)

        # Powerups for static vasp calc
        # user vasp input settings isn't working...
        wf = add_modify_incar(wf, modify_incar_params={"incar_update":
            {"ISYM": 0, "LASPH": ".TRUE.", "ADDGRID": ".TRUE.", 
            "PREC": "Accurate", "LMAXMIX": 4}})

        # Add metadata
        wf = add_additional_fields_to_taskdocs(wf, {"wf_meta": self.wf_meta})
        formula = self.matched_structures[0].composition.reduced_formula
        wf.name = "{} - Exchange".format(formula)

        tag = "magnetic_orderings group: >>{}<<".format(self.uuid)
        wf = add_tags(wf, [tag, self.ordered_structure_origins])

        return wf
