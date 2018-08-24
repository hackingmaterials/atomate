# coding: utf-8

from __future__ import absolute_import, division, print_function, \
    unicode_literals


# Defines standardized Fireworks that can be chained easily to perform various
# sequences of QChem calculations.


from fireworks import Firework

from atomate.qchem.firetasks.parse_outputs import *
from atomate.qchem.firetasks.run_calc import *
from atomate.qchem.firetasks.write_inputs import *

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "5/23/18"
__credits__ = "Brandon Wood, Shyam Dwaraknath"


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
            qchem_input_params (dict): Specify kwargs for instantiating the input set parameters.
                                       For example, if you want to change the DFT_rung, you should
                                       provide: {"DFT_rung": ...}. Defaults to None.
            db_file (str): Path to file specifying db credentials to place output parsing.
            parents ([Firework]): Parents of this particular Firework.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """

        qchem_input_params = qchem_input_params or {}
        t = []
        t.append(
            WriteInputFromIOSet(
                molecule=molecule,
                qchem_input_set="OptSet",
                input_file=input_file,
                qchem_input_params=qchem_input_params))
        t.append(
            RunQChemCustodian(
                qchem_cmd=qchem_cmd,
                multimode=multimode,
                input_file=input_file,
                output_file=output_file,
                max_cores=max_cores,
                job_type="normal"))
        t.append(
            QChemToDb(
                db_file=db_file,
                input_file=input_file,
                output_file=output_file,
                additional_fields={"task_label": name}))
        super(OptimizeFW, self).__init__(
            t,
            parents=parents,
            name=name,
            **kwargs)


class FrequencyFlatteningOptimizeFW(Firework):
    def __init__(self,
                 molecule=None,
                 name="frequency flattening structure optimization",
                 qchem_cmd="qchem",
                 multimode="openmp",
                 input_file="mol.qin",
                 output_file="mol.qout",
                 max_cores=32,
                 qchem_input_params=None,
                 max_iterations=10,
                 max_molecule_perturb_scale=0.3,
                 reversed_direction=False,
                 db_file=None,
                 parents=None,
                 **kwargs):
        """
        Iteratively optimize the given structure and flatten imaginary frequencies to ensure that
        the resulting structure is a true minima and not a saddle point.

        Args:
            molecule (Molecule): Input molecule.
            name (str): Name for the Firework.
            qchem_cmd (str): Command to run QChem. Defaults to qchem.
            multimode (str): Parallelization scheme, either openmp or mpi.
            input_file (str): Name of the QChem input file. Defaults to mol.qin.
            output_file (str): Name of the QChem output file. Defaults to mol.qout.
            max_cores (int): Maximum number of cores to parallelize over. Defaults to 32.
            qchem_input_params (dict): Specify kwargs for instantiating the input set parameters.
                                       For example, if you want to change the DFT_rung, you should
                                       provide: {"DFT_rung": ...}. Defaults to None.
            max_iterations (int): Number of perturbation -> optimization -> frequency
                                  iterations to perform. Defaults to 10.
            max_molecule_perturb_scale (float): The maximum scaled perturbation that can be
                                                applied to the molecule. Defaults to 0.3.
            reversed_direction (bool): Whether to reverse the direction of the vibrational
                                       frequency vectors. Defaults to False.
            db_file (str): Path to file specifying db credentials to place output parsing.
            parents ([Firework]): Parents of this particular Firework.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """

        qchem_input_params = qchem_input_params or {}
        t = []
        t.append(
            WriteInputFromIOSet(
                molecule=molecule,
                qchem_input_set="OptSet",
                input_file=input_file,
                qchem_input_params=qchem_input_params))
        t.append(
            RunQChemCustodian(
                qchem_cmd=qchem_cmd,
                multimode=multimode,
                input_file=input_file,
                output_file=output_file,
                max_cores=max_cores,
                job_type="opt_with_frequency_flattener",
                max_iterations=max_iterations,
                max_molecule_perturb_scale=max_molecule_perturb_scale,
                reversed_direction=reversed_direction))
        t.append(
            QChemToDb(
                db_file=db_file,
                input_file=input_file,
                output_file=output_file,
                additional_fields={
                    "task_label": name,
                    "special_run_type": "frequency_flattener"
                }))
        super(FrequencyFlatteningOptimizeFW, self).__init__(
            t,
            parents=parents,
            name=name,
            **kwargs)
