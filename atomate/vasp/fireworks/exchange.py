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
        c=None,
        structures=None,
        energies=None,
    ):
        """
        Takes a set of low-energy magnetic orderings and energies and maps
        them to a Heisenberg Model to compute exchange params.

        * heisenberg_settings: 
            cutoff (float): Starting point for nearest neighbor search.
            tol (float): Tolerance for equivalent NN bonds.
            average (bool): Compute only <J>.

        Args:
            wf_uuid (int): Unique id for record keeping.
            parent_structure (Structure): Magnetic ground state.
            parents (FireWorks): Parent FWs.
            db_file (str): Path to file containing db credentials.
            heisenberg_settings (dict): A config dict for Heisenberg model 
                mapping, detailed above.
            name (str): Labels the FW.
            c (dict): Config dict.
            structures (list): Magnetic structures.
            energies (list): Total energies of magnetic structures.

        """

        cutoff = heisenberg_settings["cutoff"]
        tol = heisenberg_settings["tol"]
        average = heisenberg_settings["average"]

        fw_name = "%s %s" % (parent_structure.composition.reduced_formula, name)

        additional_fields = {
            "task_label": fw_name,
            "exchange": {
                "calc_type": name,
                "wf_uuids": [],
                "_source_wf_uuid": wf_uuid,
            },
        }

        tasks = []

        # Generate a HeisenbergModel with only <J> exchange
        if average:
            tasks.append(
                HeisenbergModelMapping(
                    db_file=db_file,
                    wf_uuid=wf_uuid,
                    parent_structure=parent_structure,
                    cutoff=cutoff,
                    tol=tol,
                    average=average,
                    structures=structures,
                    energies=energies,
                ))
            tasks.append(
                HeisenbergModelToDb(
                    db_file=db_file,
                    wf_uuid=wf_uuid,
                    parent_structure=parent_structure,
                    cutoff=cutoff,
                    tol=tol,
                    average=average)
            )

        else:
            end_cutoff = 5.0  # Larger than this is unreasonable

            for coff in np.linspace(cutoff, end_cutoff, 5):
                tasks.append(
                    HeisenbergModelMapping(
                        db_file=db_file,
                        wf_uuid=wf_uuid,
                        parent_structure=parent_structure,
                        cutoff=coff,
                        tol=tol,
                        average=average,
                    ))
                tasks.append(
                    HeisenbergModelToDb(
                        db_file=db_file,
                        wf_uuid=wf_uuid,
                        parent_structure=parent_structure,
                        cutoff=cutoff,
                        tol=tol,
                        average=average)
                )

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
        average=True,
    ):
        """Run Vampire Monte Carlo from a HeisenbergModel.

        * mc_settings:
            mc_box_size (float): MC simulation box size in nm.
            equil_timesteps (int): Number of MC equilibration moves.
            mc_timesteps (int): Number of MC moves for averaging.

        Args:
            wf_uuid (int): Unique id for record keeping.
            parent_structure (Structure): Magnetic ground state.
            parents (FireWorks): Parent FWs.
            db_file (str): Path to file containing db credentials.
            mc_settings (dict): A configuration dict for monte carlo.
            name (str): Labels the FW.
            average (bool): Use only <J> exchange param.

        """

        fw_name = "%s %s" % (parent_structure.composition.reduced_formula, name)

        additional_fields = {
            "task_label": fw_name,
            "exchange": {
                "calc_type": name,
                "wf_uuids": [],
                "_source_wf_uuid": wf_uuid,
            },
        }

        tasks = []
        tasks.append(
            HeisenbergConvergence(
                db_file=db_file,
                wf_uuid=wf_uuid,
                parent_structure=parent_structure,
                average=average,
            )
        )

        tasks.append(
            VampireMC(
                db_file=db_file,
                wf_uuid=wf_uuid,
                parent_structure=parent_structure,
                mc_settings=mc_settings,
                average=average,
            ))
        tasks.append(
            VampireToDb(
                db_file=db_file,
                wf_uuid=wf_uuid,
                parent_structure=parent_structure,
            ),

        )

        super().__init__(tasks=tasks, name=fw_name, parents=parents)
