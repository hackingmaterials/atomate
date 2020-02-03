# coding: utf-8

from fireworks import Workflow
from atomate.qchem.fireworks.core import FrequencyFlatteningOptimizeFW, CubeAndCritic2FW
from atomate.utils.utils import get_logger

__author__ = "Samuel Blau"

logger = get_logger(__name__)


def get_wf_FFopt_and_critic(molecule,
                            suffix,
                            qchem_input_params=None,
                            db_file=">>db_file<<",
                            **kwargs):
    """
    """

    # FFopt
    fw1 = FrequencyFlatteningOptimizeFW(
         molecule=molecule,
         name="{}:{}".format(molecule.composition.alphabetical_formula, "FFopt_" + suffix),
         qchem_cmd=">>qchem_cmd<<",
         max_cores=">>max_cores<<",
         qchem_input_params=qchem_input_params,
         linked=True,
         db_file=">>db_file<<"
    )

    # Critic
    fw2 = CubeAndCritic2FW(
         name="{}:{}".format(molecule.composition.alphabetical_formula, "CC2_" + suffix),
         qchem_cmd=">>qchem_cmd<<",
         max_cores=">>max_cores<<",
         qchem_input_params=qchem_input_params,
         db_file=">>db_file<<",
         parents=fw1)
    fws = [fw1, fw2]

    wfname = "{}:{}".format(molecule.composition.alphabetical_formula, "FFopt_CC2_WF_" + suffix)

    return Workflow(fws, name=wfname, **kwargs)
