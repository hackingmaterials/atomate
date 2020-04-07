from fireworks import Firework

from pymatgen.analysis.magnetism.heisenberg import HeisenbergMapper

from atomate.vasp.firetasks.parse_outputs import JsonToDb

from atomate.vasp.firetasks.exchange import (
    HeisenbergModelMapping,
    HeisenbergModelToDb,
    HeisenbergConvergence,
    VampireMC,
    VampireToDb,
)

from atomate.vasp.config import VASP_CMD, DB_FILE

import numpy as np

__author__ = "Nathan C. Frey"
__email__ = "ncfrey@lbl.gov"


class HeisenbergModelFW(Firework):
    def __init__(
        self,
        wf_uuid,
        parent_structure,
        parents,
        db_file=DB_FILE,
        heisenberg_settings=None,
        name="heisenberg model",
        structures=None,
        energies=None,
    ):
        """
        Takes a set of low-energy magnetic orderings and energies and maps
        them to a Heisenberg Model to compute exchange params.

        * heisenberg_settings: 
            cutoff (float): Starting point for nearest neighbor search.
            tol (float): Tolerance for equivalent NN bonds.

        Args:
            wf_uuid (int): Unique id for record keeping.
            parent_structure (Structure): Magnetic ground state.
            parents (FireWorks): Parent FWs.
            db_file (str): Path to file containing db credentials.
            heisenberg_settings (dict): A config dict for Heisenberg model 
                mapping, detailed above.
            name (str): Labels the FW.
            structures (list): Magnetic structures.
            energies (list): Total energies of magnetic structures.

        TODO:
            * Test a range of nn cutoffs and add convergence check.

        """

        heisenberg_settings = heisenberg_settings or {}

        fw_name = "%s %s" % (parent_structure.composition.reduced_formula, name)

        additional_fields = {
            "task_label": fw_name,
            "exchange": {"calc_type": name, "wf_uuids": [], "_source_wf_uuid": wf_uuid},
        }

        tasks = []

        tasks.append(
            HeisenbergModelMapping(
                structures=structures,
                energies=energies,
                heisenberg_settings=heisenberg_settings,
            )
        )
        tasks.append(HeisenbergModelToDb(db_file=db_file, wf_uuid=wf_uuid))

        super().__init__(tasks=tasks, name=fw_name, parents=parents)


class VampireCallerFW(Firework):
    def __init__(
        self,
        wf_uuid,
        parent_structure,
        parents,
        db_file=DB_FILE,
        mc_settings=None,
        name="vampire caller",
    ):
        """Run Vampire Monte Carlo from a HeisenbergModel.

        * mc_settings:
            mc_box_size (float): MC simulation box size in nm.
            equil_timesteps (int): Number of MC equilibration moves.
            mc_timesteps (int): Number of MC moves for averaging.
            avg (bool): Compute only <J>.

        Args:
            wf_uuid (int): Unique id for record keeping.
            parent_structure (Structure): Magnetic ground state.
            parents (FireWorks): Parent FWs.
            db_file (str): Path to file containing db credentials.
            mc_settings (dict): A configuration dict for monte carlo.
            name (str): Labels the FW.

        """

        fw_name = "%s %s" % (parent_structure.composition.reduced_formula, name)

        additional_fields = {
            "task_label": fw_name,
            "exchange": {"calc_type": name, "wf_uuids": [], "_source_wf_uuid": wf_uuid},
        }

        tasks = []
        # tasks.append(
        #     HeisenbergConvergence(
        #         db_file=db_file,
        #         wf_uuid=wf_uuid,
        #         average=average,
        #     )
        # )

        tasks.append(
            VampireMC(db_file=db_file, wf_uuid=wf_uuid, mc_settings=mc_settings)
        )
        tasks.append(VampireToDb(db_file=db_file, wf_uuid=wf_uuid))

        super().__init__(tasks=tasks, name=fw_name, parents=parents)
