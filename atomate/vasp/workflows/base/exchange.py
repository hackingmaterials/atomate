# coding: utf-8

import os
import numpy as np

from atomate.vasp.fireworks.core import StaticFW
from fireworks import Workflow, Firework
from atomate.vasp.powerups import (
    add_tags,
    add_additional_fields_to_taskdocs,
    add_wf_metadata,
    add_common_powerups,
)
from atomate.vasp.workflows.base.core import get_wf
from atomate.vasp.firetasks.parse_outputs import (
    MagneticDeformationToDB,
    MagneticOrderingsToDB,
)

from atomate.vasp.fireworks.exchange import HeisenbergModelFW

from uuid import uuid4
from copy import deepcopy
from atomate.utils.utils import get_logger

logger = get_logger(__name__)

from atomate.vasp.config import VASP_CMD, DB_FILE, ADD_WF_METADATA
from atomate.vasp.workflows.presets.scan import wf_scan_opt

from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.analysis.structure_matcher import StructureMatcher, ElementComparator
from pymatgen.analysis.magnetism.analyzer import CollinearMagneticStructureAnalyzer, MagneticStructureEnumerator

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
        name='Exchange WF'
        ):
        """Workflow for computing exchange parameters.

        This workflow takes a set of magnetic orderings and their energies from MagneticOrderingsWF and runs static calculations of magnetic orderings in the ground state geometry. The static orderings/energies are then fit to a classical Heisenberg Hamiltonian to compute exchange parameters. The critical temperature can then be calculated with Monte Carlo.

        Args:
            magnetic_structures (list): Structure objects with the 'magmom'
        site property.
            energies (list): Total energies in eV.
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

        self.matched_structures = self.get_commensurate_orderings(magnetic_structures, energies)

    @staticmethod
    def get_commensurate_orderings(magnetic_structures, energies):
        """Generate supercells for static calculations.

        From the ground state magnetic ordering and first excited state, generate supercells with magnetic orderings.

        If the gs is FM, match supercells to the 1st es.
        If the gs is AFM, FM is included by default, 1st es may not be.
        If the gs is FiM, 1st es may not be captured.

        Args:
            magnetic_structures (list): Structures.
            energies (list): Energies.

        Returns:
            matched_structures (list): Commensurate supercells for static 
        calculations.

        """

        # Sort by energies
        ordered_structures = [
            s for _, s in sorted(zip(energies, magnetic_structures), reverse=False)
        ]
        ordered_energies = sorted(energies, reverse=False)

        # Ground state and 1st excited state
        gs_struct = ordered_structures[0]
        es_struct = ordered_structures[1]

        cmsa = CollinearMagneticStructureAnalyzer(es_struct, threshold=0.0, make_primitive=False)
        es_ordering = cmsa.ordering.value

        cmsa = CollinearMagneticStructureAnalyzer(gs_struct, threshold=0.0, make_primitive=False)
        gs_ordering = cmsa.ordering.value
        
        # FM gs will always be commensurate so we match to the 1st es
        if gs_ordering == 'FM':
            enum_struct = es_struct
        elif gs_ordering in ['FiM', 'AFM']:
            enum_struct = gs_struct
            
        mse = MagneticStructureEnumerator(enum_struct, transformation_kwargs={'min_cell_size': 1, 'max_cell_size': 2, 'check_ordered_symmetry': False})

        # Get commensurate supercells
        cmsa = CollinearMagneticStructureAnalyzer(enum_struct, threshold=0.0, make_primitive=False)
        enum_index = [cmsa.matches_ordering(s) for s in mse.ordered_structures].index(True)
        matched_structures = [mse.ordered_structures[enum_index]]

        sm = StructureMatcher(primitive_cell=False, attempt_supercell=True, 
                        comparator=ElementComparator())

        for s in mse.ordered_structures:
            s2 = sm.get_s2_like_s1(enum_struct, s) 
            if s2 is not None:
                matched_structures.append(s2)

        return matched_structures           


    def get_wf(self, num_orderings_hard_limit=16, c=None):
        """Retrieve Fireworks workflow.

        Args:
            num_orderings_hard_limit: will make sure total number of magnetic
        orderings does not exceed this number even if there are extra orderings
        of equivalent symmetry
            c: additional config dict (as used elsewhere in atomate)

        Returns: 
            wf (Workflow): Static magnetic orderings + Heisenberg Model + Vampire Monte Carlo.

        TODO:
            * Add static SCAN option
            * Make FWs modular so the wflow can be started from any step
            * Allow for user inputs to HeisenbergMapper and VampireCaller
            * Update VampireCaller to take hmodel INSTEAD of structures and energies
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

            cmsa = CollinearMagneticStructureAnalyzer(struct, threshold=0.0, make_primitive=False)

            name = " ordering {} {} -".format(idx, cmsa.ordering.value)

            static_fws.append(StaticFW(struct, vasp_cmd=c["VASP_CMD"], db_file=c["DB_FILE"], name=name + '_static'))

            analysis_parents.append(static_fws[-1])

        # Magnetic ordering analysis that generates 'magnetic_orderings'
        # collection
        fw_analysis = Firework(
            MagneticOrderingsToDB(
                db_file=c["DB_FILE"],
                wf_uuid=self.uuid,
                name="MagneticOrderingsToDB",
                parent_structure=self.matched_structures[0],
                perform_bader=False,
                scan=False,
            ),
            name="Magnetic Exchange Analysis",
            parents=analysis_parents,
            spec={"_allow_fizzled_parents": True},
        )

        # Heisenberg model mapping
        heisenberg_model_fw = HeisenbergModelFW(exchange_wf_uuid=self.uuid, parent_structure=self.matched_structures[0], parents=[fw_analysis], db_file=c['DB_FILE'])

        # Vampire Monte Carlo
        vampire_fw = VampireCallerFW(exchange_wf_uuid=self.uuid, parent_structure=self.matched_structures[0], parents=[heisenberg_model_fw], db_file=c['DB_FILE'])

        # Chain FWs together
        fws = static_fws + [fw_analysis, heisenberg_model_fw, vampire_fw]
        wf = Workflow(fws)
        wf.name = "Exchange"

        return wf



        








