# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the workflow for computing the Raman spectra.
"""

from fireworks import Firework, Workflow

from matmethods.utils.utils import get_logger
from matmethods.common.firetasks.glue_tasks import PassCalcLocs
from matmethods.vasp.fireworks.core import OptimizeFW, LepsFW, RamanFW
from matmethods.vasp.firetasks.glue_tasks import PassNormalmodesTask, CopyVaspOutputs
from matmethods.vasp.firetasks.parse_outputs import RamanSusceptibilityTensorToDbTask

from pymatgen.io.vasp.sets import MPRelaxSet

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


def get_wf_raman_spectra(structure, modes=(0, 1), step_size=0.01, vasp_cmd="vasp", db_file=None):
    """
    Raman spectra workflow:
        Calculation of phonon normal modes followed by computation of dielectric constant for
        structures displaced along the normal modes. Finally the dieledctric constants for each
        displacement is used to compute the Raman susceptibility tensor using finite difference(
        central difference scheme).

    Args:
        structure (Structure): Input structure
        vasp_input_set (DictVaspInputSet): Vasp input set for the first firework(structure optimization)
        modes (tuple/list): list of modes for which the raman spectra need to be calculated.
        step_size (float): site displacement along the normal mode in Angstroms
        vasp_cmd (str): command to run
        db_file (str): path to file containing the database credentials.

    Returns:
        Workflow
    """
    vis = MPRelaxSet(structure, force_gamma=True)
    # displacements in + and - direction along the normal mode so that the central difference scheme
    # can be used for the evaluation of Raman tensor (derivative of epsilon wrt displacement)
    displacements = [-step_size, step_size]

    fws = []

    # Structure optimization
    fw_opt = OptimizeFW(structure=structure, vasp_input_set=vis,  ediffg=-0.05, vasp_cmd=vasp_cmd,
                        db_file=db_file)
    fws.append(fw_opt)

    # Static run compute the normal modes
    fw_leps = LepsFW(structure=structure, vasp_cmd=vasp_cmd, db_file=db_file, parents=fw_opt)
    fws.append(fw_leps)

    # Extract and pass normal modes
    # why this firework: avoid parsing the xml in all subsequent fireworks,
    # also the normal modes might be needed in the final evaluation step
    fw_nm = Firework([CopyVaspOutputs(calc_loc=True, additional_files=["CHGCAR"], contcar_to_poscar=True),
                      PassNormalmodesTask(),
                      PassCalcLocs(name="pass normal modes")],
                     parents=fw_leps,
                     name="{}-{}".format(structure.composition.reduced_formula, "pass normal modes"))
    fws.append(fw_nm)

    # Static run to compute epsilon for each mode and displacement along that mode.
    fws_nm_disp = []
    for mode in modes:
        for disp in displacements:
            fws_nm_disp.append(RamanFW(structure, mode, disp, parents=fw_nm, vasp_cmd=vasp_cmd, db_file=db_file))
    fws.extend(fws_nm_disp)

    # Compute the Raman susceptibility tensor
    fw_analysis = Firework(RamanSusceptibilityTensorToDbTask(modes=modes, displacements=displacements),
                           parents=fws_nm_disp,
                           name="{}-{}".format(structure.composition.reduced_formula, "raman analysis"))
    fws.append(fw_analysis)

    wfname = "{}:{}".format(structure.composition.reduced_formula, "raman spectra")
    return Workflow(fws, name=wfname)


if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest

    structure = PymatgenTest.get_structure("Si")
    wf = get_wf_raman_spectra(structure)
