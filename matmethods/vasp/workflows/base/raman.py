# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the  workflow for Raman spectra
"""

import numpy as np
from numpy.linalg import norm

from fireworks import FireTaskBase, Firework, FWAction, Workflow
from fireworks.utilities.fw_utilities import explicit_serialize

from matmethods.utils.utils import get_logger
from matmethods.vasp.fireworks.core import OptimizeFW, LepsFW, TransmuterFW
from matmethods.vasp.firetasks.glue_tasks import CopyVaspOutputs
from matmethods.common.firetasks.glue_tasks import PassCalcLocs
from matmethods.vasp.firetasks.parse_outputs import VaspToDbTask
from matmethods.vasp.firetasks.run_calc import RunVaspCustodian
from matmethods.vasp.firetasks.write_inputs import *

from pymatgen.transformations.site_transformations import TranslateSitesTransformation
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.sets import MPRelaxSet

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


@explicit_serialize
class WriteNormalmodeTranslationTask(FireTaskBase):
    """
    """

    required_params = ["mode", "displacement", "vasp_input_set"]

    def run_task(self, fw_spec):
        vrun = Vasprun('vasprun.xml.gz')
        normalmode_eigenvecs = vrun.normalmode_eigenvecs
        nmodes, natoms, _ = normalmode_eigenvecs.shape
        # normalize the eigen vectors
        for i in range(nmodes):
            for j in range(natoms):
                normalmode_eigenvecs[i, j, :] = normalmode_eigenvecs[i, j, :] / norm(
                    normalmode_eigenvecs[i, j, :])
        structure = vrun.final_structure.copy()
        vis_cls = load_class("pymatgen.io.vasp.sets", self["vasp_input_set"])
        # add transmuter firework for the given mode and displacement
        normalmode_translation = structure.frac_coords + normalmode_eigenvecs[self["mode"], :, :] * self["displacement"]
        transformation = TranslateSitesTransformation(indices_to_move=range(len(structure)),
                                                    translation_vector=normalmode_translation.tolist(),
                                                    vector_in_frac_coords=True)
        ts = TransformedStructure(structure)
        transmuter = StandardTransmuter([ts], [transformation])
        vis = vis_cls(transmuter.transformed_structures[-1].final_structure, **self.get("vasp_input_params", {}))
        vis.write_input(".")


@explicit_serialize
class PassEpsilonTask(FireTaskBase):
    """
    pass mode, displacement and the corresponding epsilon
    """

    required_params = ["mode", "displacement"]

    def run_task(self, fw_spec):
        vrun = Vasprun('vasprun.xml.gz')
        epsilon_static = vrun.epsilon_static
        return FWAction(mod_spec=[{
            '_set': {'raman_{}_{}'.format(self["mode"], self["displacement"]): epsilon_static}
        }])


@explicit_serialize
class RamanAnalysisTask(FireTaskBase):
    """
    derivative of epsilon_static wrt position along the normal mode
    """

    required_params = ["modes", "displacements"]

    def run_task(self, fw_spec):
        pass


class RamanFW(Firework):
    def __init__(self, structure, mode, displacement, name="raman spectra", vasp_cmd="vasp",
                 db_file=None, parents=None, **kwargs):
        """
        Firework for raman spectra calculation.

        Args:

        """
        t = []

        t.append(CopyVaspOutputs(calc_loc=True, contcar_to_poscar=True))
        t.append(WriteNormalmodeTranslationTask(mode=mode, displacement=displacement))
        t.append(RunVaspCustodian(vasp_cmd=vasp_cmd))
        t.append(PassCalcLocs(name=name))
        t.append(PassEpsilonTask(mode=mode, displacement=displacement))
        t.append(VaspToDbTask(db_file=db_file,
                              additional_fields={"task_label": name,
                                                 "transmuter": {"transformations": transformations,
                                                                "transformation_params": transformation_params}
                                                 }))
        super(RamanFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)


def get_wf_raman_spectra(structure, vasp_input_set=None, modes=(0,1), step_size=0.01, vasp_cmd="vasp",
                         db_file=None):
    """
    Raman spectra workflow

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
    fw_leps = LepsFW(structure=structure, vasp_cmd=vasp_cmd, db_file=db_file, parents=[fw_opt])
    fws.append(fw_leps)

    # add raman firework for each mode and displacement along that mode
    for mode in modes:
        for disp in displacements:
            fws.append(RamanFW(structure, mode, disp, vasp_cmd=vasp_cmd, db_file=db_file, parents=[fw_leps]))

    # analysis: compute the intensity
    #fws.append(Firework(RamanAnalysisTask(modes=modes, displacements=displacements)))

    wfname = "{}:{}".format(structure.composition.reduced_formula, "raman spectra")

    return Workflow(fws, name=wfname)


if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest

    structure = PymatgenTest.get_structure("Si")
    wf = get_wf_raman_spectra(structure)
