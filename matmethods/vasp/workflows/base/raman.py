# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the workflow for computing the Raman spectra.
"""

from collections import defaultdict

import numpy as np
from numpy.linalg import norm

from fireworks import FireTaskBase, Firework, FWAction, Workflow
from fireworks.utilities.fw_utilities import explicit_serialize

from matmethods.utils.utils import get_logger
from matmethods.vasp.fireworks.core import OptimizeFW, LepsFW
from matmethods.vasp.firetasks.glue_tasks import CopyVaspOutputs
from matmethods.common.firetasks.glue_tasks import PassCalcLocs
from matmethods.vasp.firetasks.parse_outputs import VaspToDbTask
from matmethods.vasp.firetasks.run_calc import RunVaspCustodian

from pymatgen.alchemy.materials import TransformedStructure
from pymatgen.alchemy.transmuters import StandardTransmuter
from pymatgen.transformations.site_transformations import TranslateSitesTransformation
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


@explicit_serialize
class WriteNormalmodeDisplacementTask(FireTaskBase):
    """
    write the vasp input set corresponding to the structure displaced along the normal mode.
    """

    required_params = ["mode", "displacement", "vasp_input_set"]
    optional_params = ["vasp_input_params"]

    def run_task(self, fw_spec):
        vrun = Vasprun('vasprun.xml.gz')
        structure = vrun.final_structure.copy()
        normalmode_eigenvecs = vrun.normalmode_eigenvecs
        nmodes, natoms, _ = normalmode_eigenvecs.shape
        # normalize the eigen vectors
        for i in range(nmodes):
            for j in range(natoms):
                normalmode_eigenvecs[i, j, :] = normalmode_eigenvecs[i, j, :] / norm(normalmode_eigenvecs[i, j, :])

        # displace the sites along the given normal mode
        normalmode_translation = structure.cart_coords + normalmode_eigenvecs[self["mode"], :, :] * self["displacement"]
        transformation = TranslateSitesTransformation(range(len(structure)), normalmode_translation,
                                                      vector_in_frac_coords=True)
        ts = TransformedStructure(structure)
        transmuter = StandardTransmuter([ts], [transformation])

        # write the static vasp input set corresponding to the transmuted structure
        vis = self["vasp_input_set"].__class__(transmuter.transformed_structures[-1].final_structure,
                                               **self.get("vasp_input_params", {}))
        vis.write_input(".")


@explicit_serialize
class PassEpsilonTask(FireTaskBase):
    """
    pass the epsilon corresponding to the mode and displacement.
    """

    required_params = ["mode", "displacement"]

    def run_task(self, fw_spec):
        vrun = Vasprun('vasprun.xml.gz')
        epsilon_static = vrun.epsilon_static
        epsilon_dict = {"mode": self["mode"],
                        "displacement": self["displacement"],
                        "epsilon": epsilon_static}
        return FWAction(mod_spec=[{
            '_set': {
                'raman_epsilon->{}_{}'.format(str(self["mode"]), str(self["displacement"])): epsilon_dict
            }
        }])


@explicit_serialize
class RamanSusceptibilityTensorTask(FireTaskBase):
    """
    finite difference derivative of epsilon_static wrt position along the normal mode
    --> raman susceptibilty tensor for each mode. See: 10.1103/PhysRevB.63.094305

    """

    required_params = ["modes", "displacements"]

    def run_task(self, fw_spec):
        mode_disps = fw_spec["raman_epsilon"].keys()
        # store the dispalcement & epsilon for each mode in a dictionary
        modes_eps_dict = defaultdict(list)
        for md in mode_disps:
            modes_eps_dict[str(fw_spec["raman_epsilon"][md]["mode"])].append(
                [fw_spec["raman_epsilon"][md]["displacement"],
                 fw_spec["raman_epsilon"][md]["epsilon"]])

        # raman tensor = finite difference derivative of epsilon wrt displacement.
        for k, v in modes_eps_dict.items():
            modes_eps_dict["raman"] = (np.array(v[0][1]) - np.array(v[1][1]))/(v[0][0] - v[1][0])


class RamanFW(Firework):
    def __init__(self, structure, mode, displacement, parents, vasp_input_set=None, name="raman spectra",
                 vasp_cmd="vasp", db_file=None, **kwargs):
        """
        Firework for raman spectra calculation.

        Args:
            structure (Structure): Input structure.
            mode (int): index of the normal mode
            displacement (float): displacement along the normal mode in Angstroms
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            vasp_input_set VaspInputSet): input set to use (for jobs w/no parents)
                Defaults to MPStaticSet() if None.
            name (str): Name for the Firework.
            vasp_cmd (str): Command to run vasp.
            db_file (str): Path to file specifying db credentials.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.

        """
        t = []
        vasp_input_set = vasp_input_set or MPStaticSet(structure)
        t.append(CopyVaspOutputs(calc_loc=True, contcar_to_poscar=True))
        t.append(WriteNormalmodeDisplacementTask(mode=mode, displacement=displacement, vasp_input_set=vasp_input_set))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd))
        t.append(PassCalcLocs(name=name))
        t.append(PassEpsilonTask(mode=mode, displacement=displacement))
        t.append(VaspToDbTask(db_file=db_file, additional_fields={"task_label": name}))
        super(RamanFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)


def get_wf_raman_spectra(structure, vasp_input_set=None, modes=(0, 1), step_size=0.01, vasp_cmd="vasp",
                         db_file=None):
    """
    Raman spectra workflow:
        Calculation of phonon normal modes followed by computation of dielectric constant for
        structures displaced along the normal modes. Finally the dieledctric constants for each
        displacement is used to compute the Raman susceptibility tensor using finite difference(
        central difference scheme).

    Returns:
        Workflow
    """
    vis = vasp_input_set or MPRelaxSet(structure, force_gamma=True)
    # displacements in + and - direction along the normal mode so that the central difference scheme
    # can be used
    displacements = [-step_size, step_size]

    fws = []

    # optimize
    fw_opt = OptimizeFW(structure=structure, vasp_input_set=vis, vasp_cmd=vasp_cmd, db_file=db_file)
    fws.append(fw_opt)

    # leps firework to compute the normal modes
    fw_leps = LepsFW(structure=structure, vasp_cmd=vasp_cmd, db_file=db_file, parents=fw_opt)
    fws.append(fw_leps)

    # add raman firework for each mode and displacement along that mode
    for mode in modes:
        for disp in displacements:
            fws.append(RamanFW(structure, mode, disp, fw_leps, vasp_cmd=vasp_cmd, db_file=db_file))

    # analysis: compute the intensity
    #fws.append(Firework(RamanSusceptibilityTensorTask(modes=modes, displacements=displacements)))

    wfname = "{}:{}".format(structure.composition.reduced_formula, "raman spectra")
    return Workflow(fws, name=wfname)


if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest

    structure = PymatgenTest.get_structure("Si")
    wf = get_wf_raman_spectra(structure)
