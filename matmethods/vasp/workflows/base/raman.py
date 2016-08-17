# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines the  workflow for raman spectra
"""

from fireworks import FireTaskBase, Firework, FWAction, Workflow
from fireworks.utilities.fw_utilities import explicit_serialize

from matmethods.utils.utils import get_logger
from matmethods.vasp.fireworks.core import OptimizeFW, LepsFW, TransmuterFW

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.sets import MPRelaxSet

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

logger = get_logger(__name__)


@explicit_serialize
class NormalmodeDeformationsTask(FireTaskBase):
    """
    """
    required_params = ["modes", "displacements"]
    def run_task(self, fw_spec):
        vrun = Vasprun('vasprun.xml.gz')
        structure = vrun.final_structure
        normalmode_eigenvecs = vrun.normalmode_eigenvecs
        fws = []
        # add transmuter firework for each mode and displacement
        for idx in self["modes"]:
            for disp in self["displacements"]:
                normalmode_defo = None # define the eigen vec based deformation here, normalmode_eigenvecs[idx]
                fw = TransmuterFW(name="normal mode deformation",
                                 structure=structure,
                                 transformations=['DeformStructureTransformation'],
                                 transformation_params=[{"deformation": normalmode_defo.tolist()}],
                                 db_file=self["db_file"],
                                 vasp_cmd=self["vasp_cmd"],
                                 )
                fws.append(fw)
        return FWAction(additions=fws)


@explicit_serialize
class RamanAnalysisTask(FireTaskBase):
    """
    derivative of epsilon_static wrt position along the normal mode
    """
    required_params = ["modes", "displacements"]
    def run_task(self, fw_spec):
        pass


def get_wf_raman_spectra(structure, vasp_input_set=None, modes=(0,1), step_size=0.01,
                         vasp_cmd="vasp", db_file=None):
    """

    Returns:
        Workflow
    """
    vis = vasp_input_set or MPRelaxSet(structure, force_gamma=True)
    # displacements in + and - direction along the normal mode so that the central difference scheme
    # can be used
    displacements = [-step_size, step_size]

    fws=[]

    # optimize
    fws.append(OptimizeFW(structure=structure, eigenvasp_input_set=vis, vasp_cmd=vasp_cmd, db_file=db_file))

    # leps + task that generates raman fireworks
    fw_leps = LepsFW(structure=structure, vasp_cmd=vasp_cmd, db_file=db_file)
    fw_leps.spec['_tasks'].append(NormalmodeDeformationsTask(modes=modes, displacements=displacements).to_dict())
    fws.append(fw_leps)

    # analysis: compute the intensity
    fws.append(Firework(RamanAnalysisTask(modes=modes, displacements=displacements)))

    wfname = "{}:{}".format(structure.composition.reduced_formula, "raman spectra")

    return Workflow(fws, name=wfname)


if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest

    structure = PymatgenTest.get_structure("Si")
    wf = get_wf_raman_spectra(structure)
