# coding: utf-8

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import warnings

"""
Defines standardized Fireworks that can be chained easily to perform various
sequences of VASP calculations.
"""

from fireworks import Firework

from pymatgen import Structure
from pymatgen.io.qchem_io.sets import *

from atomate.qchem.firetasks.parse_outputs import *
from atomate.qchem.firetasks.run_calc import *
from atomate.qchem.firetasks.write_inputs import *


class OptimizeFW(Firework):

    def __init__(self, 
                 molecule=None,
                 name="structure optimization",
                 qchem_cmd="qchem",
                 multimode="openmp",
                 input_file="mol.qin",
                 output_file="mol.qout",
                 max_cores=32,
                 qchem_input_params=None,
                 db_file=None,
                 parents=None,
                 **kwargs):
        """
        Optimize the given structure.

        Args:
            molecule (Molecule): Input molecule.
            name (str): Name for the Firework.
            qchem_cmd (str): Command to run QChem. Defaults to qchem.
            multimode (str): Parallelization scheme, either openmp or mpi.
            input_file (str): Name of the QChem input file. Defaults to mol.qin.
            output_file (str): Name of the QChem output file. Defaults to mol.qout.
            max_cores (int): Maximum number of cores to parallelize over. Defaults to 32.
            save_name (str): Name of the saved scratch directory. Defaults to "default_save_name".
            qchem_input_params (dict): Specify kwargs for instantiating the input set parameters. 
                                       For example, if you want to change the DFT_rung, you should 
                                       provide: {"DFT_rung": ...}. Defaults to None.
            db_file (str): Path to file specifying db credentials to place output parsing.
            job_type (str): custodian job type
            parents ([Firework]): Parents of this particular Firework.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """

        t = []
        t.append(WriteInputFromIOSet(molecule=molecule,
                                     qchem_input_set="OptSet",
                                     input_file=input_file,
                                     qchem_input_params=qchem_input_params))
        t.append(RunQChemCustodian(qchem_cmd=qchem_cmd,
                                   multimode=multimode,
                                   input_file=input_file,
                                   output_file=output_file,
                                   max_cores=max_cores,
                                   job_type="normal"))
        t.append(QChemToDb(db_file=db_file,
                           input_file=input_file,
                           output_file=output_file,
                           additional_fields={"task_label": name}))
        super(OptimizeFW, self).__init__(t, parents=parents,
                                         name="{}-{}".format(structure.composition.reduced_formula, name),
                                         **kwargs)

