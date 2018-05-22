# coding: utf-8

from __future__ import absolute_import, division, print_function, \
    unicode_literals

"""
This module defines the torsion potential workflow
"""

import numpy as np

from fireworks import Firework, Workflow

from atomate.utils.utils import get_logger, get_fws_and_tasks

__author__ = 'Brandon Wood'
__email__ = 'b.wood@berkeley.edu'

logger = get_logger(__name__)


def get_wf_torsion_potential(molecule, name="torsion_potential"
                            db_file=None,
                            conventional=False, order=2, vasp_input_set=None,
                            analysis=True,
                            sym_reduce=False, tag='elastic',
                            copy_vasp_outputs=False, **kwargs):
    """
    Returns a workflow to the torsion potential for a molecule.

    Firework 1 : write QChem input for an optimization,
                 run Qchem,
                 parse output and insert into db,
                 pass relaxed molecule to fw_spec and on to fw2,

    Firework 2 : rotate molecule torsion to a particular angle,
                 write QChem input for an optimization,
                 run Qchem,
                 parse output and insert into db

    last Firework : add analysis code at some point

    Args:
        molecule (Molecule): input molecule to be optimized and run.
        db_file (str): path to file containing the database credentials.
        conventional (bool): flag to convert input structure to conventional structure,
            defaults to False.
        order (int): order of the tensor expansion to be determined.  Defaults to 2 and
            currently supports up to 3.
        vasp_input_set (VaspInputSet): vasp input set to be used.  Defaults to static
            set with ionic relaxation parameters set.  Take care if replacing this,
            default ensures that ionic relaxation is done and that stress is calculated
            for each vasp run.
        analysis (bool): flag to indicate whether analysis task should be added
            and stresses and strains passed to that task
        sym_reduce (bool): Whether or not to apply symmetry reductions
        tag (str):
        copy_vasp_outputs (bool): whether or not to copy previous vasp outputs.
        kwargs (keyword arguments): additional kwargs to be passed to get_wf_deformations

    Returns:
        Workflow
    """
    fws, parents = [], []

    qchem_input_set = "OptSet"

    # Optimize the starting molecule fw1
    fws.append(OptimizeFW())

    # Loop to generate all the different rotated molecule optimizations
    for angle in angles:

        fws.append(OptimizeFW())

    wfname = "{}:{}".format(Molecule.formula_name, name)

    return Workflow(fws, name=wfname)
    


