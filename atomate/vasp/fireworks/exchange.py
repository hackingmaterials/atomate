from fireworks import Firework

from pymatgen.analysis.magnetism.heisenberg import HeisenbergMapper

from atomate.vasp.firetasks.parse_outputs import JsonToDb

from atomate.vasp.firetasks.exchange_tasks import (
    HeisenbergModelMapping,
    HeisenbergConvergence,
    VampireMC,
)

from atomate.vasp.config import VASP_CMD, DB_FILE

import numpy as np

__author__ = "Nathan C. Frey"
__email__ = "ncfrey@lbl.gov"


class HeisenbergModelFW(Firework):
    def __init__(
        self,
        exchange_wf_uuid,
        parent_structure,
        parents,
        db_file=DB_FILE,
        cutoff=3.0,
        tol=0.04,
        name="heisenberg model",
        c=None,
    ):
        """
        Takes a set of low-energy magnetic orderings and energies and maps them to a Heisenberg Model to compute exchange params.

        Args:
            exchange_wf_uuid (int): Unique id for record keeping.
            parent_structure (Structure): Magnetic ground state.
            parents (FireWorks): Parent FWs.
            db_file (str): Path to file containing db credentials.
            cutoff (float): Starting point for nearest neighbor search.
            tol (float): Tolerance for equivalent NN bonds.
            name (str): Labels the FW.
            c (dict): Config dict.

        """

        fw_name = "%s %s" % (parent_structure.composition.reduced_formula, name)

        additional_fields = {
            "task_label": fw_name,
            "exchange": {
                "calc_type": name,
                "wf_uuids": [],
                "_source_wf_uuid": exchange_wf_uuid,
            },
        }

        tasks = []
        end_cutoff = cutoff + 15

        for c in np.linspace(cutoff, end_cutoff, 9):
            tasks.append(
                HeisenbergModelMapping(
                    db_file=db_file,
                    exchange_wf_uuid=exchange_wf_uuid,
                    parent_structure=parent_structure,
                    cutoff=c,
                    tol=tol,
                )
            )

        super().__init__(tasks=tasks, name=fw_name, parents=parents)


class VampireCallerFW(Firework):
    def __init__(
        self,
        exchange_wf_uuid,
        parent_structure,
        parents,
        db_file=DB_FILE,
        name="vampire caller",
        c=None,
    ):
        """Run Vampire Monte Carlo from a HeisenbergModel.

        Args:
            exchange_wf_uuid (int): Unique id for record keeping.
            parent_structure (Structure): Magnetic ground state.
            parents (FireWorks): Parent FWs.
            db_file (str): Path to file containing db credentials.
            name (str): Labels the FW.
            c (dict): Config dict.

        """

        fw_name = "%s %s" % (parent_structure.composition.reduced_formula, name)

        additional_fields = {
            "task_label": fw_name,
            "exchange": {
                "calc_type": name,
                "wf_uuids": [],
                "_source_wf_uuid": exchange_wf_uuid,
            },
        }

        tasks = []
        tasks.append(
            HeisenbergConvergence(
                db_file=db_file,
                exchange_wf_uuid=exchange_wf_uuid,
                parent_structure=parent_structure,
            )
        )

        tasks.append(
            VampireMC(
                db_file=db_file,
                exchange_wf_uuid=exchange_wf_uuid,
                parent_structure=parent_structure,
            )
        )

        super().__init__(tasks=tasks, name=fw_name, parents=parents)
