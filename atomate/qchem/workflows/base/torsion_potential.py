# coding: utf-8

from __future__ import absolute_import, division, print_function, \
    unicode_literals

"""
This module defines the torsion potential workflow
"""

import numpy as np

from fireworks import Firework, Workflow
from atomate.qchem.fireworks.core import OptimizeFW
from atomate.utils.utils import get_logger, get_fws_and_tasks
from atomate.qchem.firetasks.geo_transformations import RotateTorsion

__author__ = 'Brandon Wood'
__email__ = 'b.wood@berkeley.edu'

logger = get_logger(__name__)


def get_wf_torsion_potential(molecule,
                             atom_indexes,
                             angles,
                             name="torsion_potential",
                             qchem_cmd="qchem",
                             multimode="openmp",
                             input_file="mol.qin",
                             output_file="mol.qout",
                             max_cores=32,
                             qchem_input_params=None,
                             db_file=None,
                             **kwargs):

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
            molecule (Molecule): Input molecule (needs to be a pymatgen molecule object)
            atom_indexes (list of ints): list of atom indexes in the torsion angle to be rotated (i.e. [6, 8, 9, 10])
            angles (list of floats): list of all the torsion angles to run
            name (str): Name for the workflow.
            qchem_cmd (str): Command to run QChem. Defaults to qchem.
            multimode (str): Parallelization scheme, either openmp or mpi.
            input_file (str): Name of the QChem input file. Defaults to mol.qin.
            output_file (str): Name of the QChem output file. Defaults to mol.qout.
            max_cores (int): Maximum number of cores to parallelize over. Defaults to 32.
            qchem_input_params (dict): Specify kwargs for instantiating the input set parameters.
                                       For example, if you want to change the DFT_rung, you should
                                       provide: {"DFT_rung": ...}. Defaults to None.
            db_file (str): Path to file specifying db credentials to place output parsing.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.

    Returns: Workflow
    """
    fws = []

    # Optimize the starting molecule fw1
    fw1 = OptimizeFW(molecule=molecule, name=name, qchem_cmd=qchem_cmd,
                     multimode=multimode, input_file=input_file, output_file=output_file,
                     max_cores=max_cores, qchem_input_params=qchem_input_params,
                     db_file=db_file, **kwargs)
    fws.append(fw1)

    # Loop to generate all the different rotated molecule optimizations
    for angle in angles:
        rot_opt_fw = OptimizeFW(name=name, qchem_cmd=qchem_cmd,
                                multimode=multimode, input_file=input_file, output_file=output_file,
                                max_cores=max_cores, qchem_input_params=qchem_input_params,
                                db_file=db_file, parents=fw1, **kwargs)
        rot_task = RotateTorsion(atom_indexes=atom_indexes, angle=angle)
        rot_opt_fw.tasks.insert(0, rot_task)
        fws.append(rot_opt_fw)

    wfname = "{}:{}".format(molecule.formula, name)

    return Workflow(fws, name=wfname)
    


