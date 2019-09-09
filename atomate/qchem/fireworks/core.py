# coding: utf-8

# Defines standardized Fireworks that can be chained easily to perform various
# sequences of QChem calculations.


from fireworks import Firework

from atomate.qchem.firetasks.parse_outputs import QChemToDb
from atomate.qchem.firetasks.run_calc import RunQChemCustodian
from atomate.qchem.firetasks.write_inputs import WriteInputFromIOSet
from atomate.qchem.firetasks.fragmenter import FragmentMolecule

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "5/23/18"
__credits__ = "Brandon Wood, Shyam Dwaraknath"


class SinglePointFW(Firework):
    def __init__(self,
                 molecule=None,
                 name="single point",
                 qchem_cmd=">>qchem_cmd<<",
                 multimode=">>multimode<<",
                 max_cores=">>max_cores<<",
                 qchem_input_params=None,
                 db_file=None,
                 parents=None,
                 **kwargs):
        """

        Args:
            molecule (Molecule): Input molecule.
            name (str): Name for the Firework.
            qchem_cmd (str): Command to run QChem. Supports env_chk.
            multimode (str): Parallelization scheme, either openmp or mpi. Supports env_chk.
            max_cores (int): Maximum number of cores to parallelize over. Supports env_chk.
            qchem_input_params (dict): Specify kwargs for instantiating the input set parameters.
                                       Basic uses would be to modify the default inputs of the set,
                                       such as dft_rung, basis_set, pcm_dielectric, scf_algorithm,
                                       or max_scf_cycles. See pymatgen/io/qchem/sets.py for default
                                       values of all input parameters. For instance, if a user wanted
                                       to use a more advanced DFT functional, include a pcm with a
                                       dielectric of 30, and use a larger basis, the user would set
                                       qchem_input_params = {"dft_rung": 5, "pcm_dielectric": 30,
                                       "basis_set": "6-311++g**"}. However, more advanced customization
                                       of the input is also possible through the overwrite_inputs key
                                       which allows the user to directly modify the rem, pcm, smd, and
                                       solvent dictionaries that QChemDictSet passes to inputs.py to
                                       print an actual input file. For instance, if a user wanted to
                                       set the sym_ignore flag in the rem section of the input file
                                       to true, then they would set qchem_input_params = {"overwrite_inputs":
                                       "rem": {"sym_ignore": "true"}}. Of course, overwrite_inputs
                                       could be used in conjuction with more typical modifications,
                                       as seen in the test_double_FF_opt workflow test.
            db_file (str): Path to file specifying db credentials to place output parsing.
            parents ([Firework]): Parents of this particular Firework.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """

        qchem_input_params = qchem_input_params or {}
        input_file="mol.qin"
        output_file="mol.qout"
        t = []
        t.append(
            WriteInputFromIOSet(
                molecule=molecule,
                qchem_input_set="SinglePointSet",
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
        super(SinglePointFW, self).__init__(
            t,
            parents=parents,
            name=name,
            **kwargs)


class OptimizeFW(Firework):
    def __init__(self,
                 molecule=None,
                 name="structure optimization",
                 qchem_cmd=">>qchem_cmd<<",
                 multimode=">>multimode<<",
                 max_cores=">>max_cores<<",
                 qchem_input_params=None,
                 db_file=None,
                 parents=None,
                 **kwargs):
        """
        Optimize the given structure.

        Args:
            molecule (Molecule): Input molecule.
            name (str): Name for the Firework.
            qchem_cmd (str): Command to run QChem. Supports env_chk.
            multimode (str): Parallelization scheme, either openmp or mpi. Defaults to openmp.
            max_cores (int): Maximum number of cores to parallelize over. Supports env_chk.
            qchem_input_params (dict): Specify kwargs for instantiating the input set parameters.
                                       Basic uses would be to modify the default inputs of the set,
                                       such as dft_rung, basis_set, pcm_dielectric, scf_algorithm,
                                       or max_scf_cycles. See pymatgen/io/qchem/sets.py for default
                                       values of all input parameters. For instance, if a user wanted
                                       to use a more advanced DFT functional, include a pcm with a
                                       dielectric of 30, and use a larger basis, the user would set
                                       qchem_input_params = {"dft_rung": 5, "pcm_dielectric": 30,
                                       "basis_set": "6-311++g**"}. However, more advanced customization
                                       of the input is also possible through the overwrite_inputs key
                                       which allows the user to directly modify the rem, pcm, smd, and
                                       solvent dictionaries that QChemDictSet passes to inputs.py to
                                       print an actual input file. For instance, if a user wanted to
                                       set the sym_ignore flag in the rem section of the input file
                                       to true, then they would set qchem_input_params = {"overwrite_inputs":
                                       "rem": {"sym_ignore": "true"}}. Of course, overwrite_inputs
                                       could be used in conjuction with more typical modifications,
                                       as seen in the test_double_FF_opt workflow test.
            db_file (str): Path to file specifying db credentials to place output parsing.
            parents ([Firework]): Parents of this particular Firework.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """

        qchem_input_params = qchem_input_params or {}
        input_file="mol.qin"
        output_file="mol.qout"
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


class FrequencyFW(Firework):
    def __init__(self,
                 molecule=None,
                 name="frequency calculation",
                 qchem_cmd=">>qchem_cmd<<",
                 multimode=">>multimode<<",
                 max_cores=">>max_cores<<",
                 qchem_input_params=None,
                 db_file=None,
                 parents=None,
                 **kwargs):
        """
        Optimize the given structure.

        Args:
            molecule (Molecule): Input molecule.
            name (str): Name for the Firework.
            qchem_cmd (str): Command to run QChem. Supports env_chk.
            multimode (str): Parallelization scheme, either openmp or mpi. Defaults to openmp.
            max_cores (int): Maximum number of cores to parallelize over. Supports env_chk.
            qchem_input_params (dict): Specify kwargs for instantiating the input set parameters.
                                       Basic uses would be to modify the default inputs of the set,
                                       such as dft_rung, basis_set, pcm_dielectric, scf_algorithm,
                                       or max_scf_cycles. See pymatgen/io/qchem/sets.py for default
                                       values of all input parameters. For instance, if a user wanted
                                       to use a more advanced DFT functional, include a pcm with a
                                       dielectric of 30, and use a larger basis, the user would set
                                       qchem_input_params = {"dft_rung": 5, "pcm_dielectric": 30,
                                       "basis_set": "6-311++g**"}. However, more advanced customization
                                       of the input is also possible through the overwrite_inputs key
                                       which allows the user to directly modify the rem, pcm, smd, and
                                       solvent dictionaries that QChemDictSet passes to inputs.py to
                                       print an actual input file. For instance, if a user wanted to
                                       set the sym_ignore flag in the rem section of the input file
                                       to true, then they would set qchem_input_params = {"overwrite_inputs":
                                       "rem": {"sym_ignore": "true"}}. Of course, overwrite_inputs
                                       could be used in conjuction with more typical modifications,
                                       as seen in the test_double_FF_opt workflow test.
            db_file (str): Path to file specifying db credentials to place output parsing.
            parents ([Firework]): Parents of this particular Firework.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """

        qchem_input_params = qchem_input_params or {}
        input_file="mol.qin"
        output_file="mol.qout"
        t = []
        t.append(
            WriteInputFromIOSet(
                molecule=molecule,
                qchem_input_set="FreqSet",
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
        super(FrequencyFW, self).__init__(
            t,
            parents=parents,
            name=name,
            **kwargs)


class FrequencyFlatteningOptimizeFW(Firework):
    def __init__(self,
                 molecule=None,
                 name="frequency flattening structure optimization",
                 qchem_cmd=">>qchem_cmd<<",
                 multimode=">>multimode<<",
                 max_cores=">>max_cores<<",
                 qchem_input_params=None,
                 max_iterations=10,
                 max_molecule_perturb_scale=0.3,
                 linked=False,
                 db_file=None,
                 parents=None,
                 **kwargs):
        """
        Iteratively optimize the given structure and flatten imaginary frequencies to ensure that
        the resulting structure is a true minima and not a saddle point.

        Args:
            molecule (Molecule): Input molecule.
            name (str): Name for the Firework.
            qchem_cmd (str): Command to run QChem. Supports env_chk.
            multimode (str): Parallelization scheme, either openmp or mpi. Supports env_chk.
            max_cores (int): Maximum number of cores to parallelize over. Supports env_chk.
            qchem_input_params (dict): Specify kwargs for instantiating the input set parameters.
                                       Basic uses would be to modify the default inputs of the set,
                                       such as dft_rung, basis_set, pcm_dielectric, scf_algorithm,
                                       or max_scf_cycles. See pymatgen/io/qchem/sets.py for default
                                       values of all input parameters. For instance, if a user wanted
                                       to use a more advanced DFT functional, include a pcm with a
                                       dielectric of 30, and use a larger basis, the user would set
                                       qchem_input_params = {"dft_rung": 5, "pcm_dielectric": 30,
                                       "basis_set": "6-311++g**"}. However, more advanced customization
                                       of the input is also possible through the overwrite_inputs key
                                       which allows the user to directly modify the rem, pcm, smd, and
                                       solvent dictionaries that QChemDictSet passes to inputs.py to
                                       print an actual input file. For instance, if a user wanted to
                                       set the sym_ignore flag in the rem section of the input file
                                       to true, then they would set qchem_input_params = {"overwrite_inputs":
                                       "rem": {"sym_ignore": "true"}}. Of course, overwrite_inputs
                                       could be used in conjuction with more typical modifications,
                                       as seen in the test_double_FF_opt workflow test.
            max_iterations (int): Number of perturbation -> optimization -> frequency
                                  iterations to perform. Defaults to 10.
            max_molecule_perturb_scale (float): The maximum scaled perturbation that can be
                                                applied to the molecule. Defaults to 0.3.
            db_file (str): Path to file specifying db credentials to place output parsing.
            parents ([Firework]): Parents of this particular Firework.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """

        qchem_input_params = qchem_input_params or {}
        input_file="mol.qin"
        output_file="mol.qout"
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
                linked=linked))
        t.append(
            QChemToDb(
                db_file=db_file,
                input_file=input_file,
                output_file=output_file,
                additional_fields={
                    "task_label": name,
                    "special_run_type": "frequency_flattener",
                    "linked": linked
                }))
        super(FrequencyFlatteningOptimizeFW, self).__init__(
            t,
            parents=parents,
            name=name,
            **kwargs)


class FragmentFW(Firework):
    def __init__(self,
                 molecule=None,
                 depth=1,
                 open_rings=True,
                 additional_charges=None,
                 do_triplets=True,
                 linked=False,
                 name="fragment and optimize",
                 qchem_input_params=None,
                 db_file=None,
                 check_db=True,
                 parents=None,
                 **kwargs):
        """
        Fragment the given structure and optimize all unique fragments

        Args:
            molecule (Molecule): Input molecule.
            depth (int): Fragmentation depth. Defaults to 1. See fragmenter firetask for more details.
            open_rings (bool): Whether or not to open any rings encountered during fragmentation.
                               Defaults to True. See fragmenter firetask for more details.
            additional_charges (list): List of additional charges besides the defaults. See fragmenter
                                       firetask for more details.
            do_triplets (bool): Whether to simulate triplets as well as singlets for molecules with an
                                even number of electrons. Defaults to True.
            name (str): Name for the Firework.
            qchem_input_params (dict): Specify kwargs for instantiating the input set parameters.
                                       Basic uses would be to modify the default inputs of the set,
                                       such as dft_rung, basis_set, pcm_dielectric, scf_algorithm,
                                       or max_scf_cycles. See pymatgen/io/qchem/sets.py for default
                                       values of all input parameters. For instance, if a user wanted
                                       to use a more advanced DFT functional, include a pcm with a
                                       dielectric of 30, and use a larger basis, the user would set
                                       qchem_input_params = {"dft_rung": 5, "pcm_dielectric": 30,
                                       "basis_set": "6-311++g**"}. However, more advanced customization
                                       of the input is also possible through the overwrite_inputs key
                                       which allows the user to directly modify the rem, pcm, smd, and
                                       solvent dictionaries that QChemDictSet passes to inputs.py to
                                       print an actual input file. For instance, if a user wanted to
                                       set the sym_ignore flag in the rem section of the input file
                                       to true, then they would set qchem_input_params = {"overwrite_inputs":
                                       "rem": {"sym_ignore": "true"}}. Of course, overwrite_inputs
                                       could be used in conjuction with more typical modifications,
                                       as seen in the test_double_FF_opt workflow test.
            db_file (str): Path to file specifying db credentials to place output parsing.
            check_db (bool): Whether or not to check the database for equivalent structures
                             before adding new fragment fireworks. Defaults to True.
            parents ([Firework]): Parents of this particular Firework.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """

        qchem_input_params = qchem_input_params or {}
        additional_charges = additional_charges or []
        t = []
        t.append(
            FragmentMolecule(
                molecule=molecule,
                depth=depth,
                open_rings=open_rings,
                additional_charges=additional_charges,
                do_triplets=do_triplets,
                linked=linked,
                qchem_input_params=qchem_input_params,
                db_file=db_file,
                check_db=check_db))
        super(FragmentFW, self).__init__(
            t,
            parents=parents,
            name=name,
            **kwargs)
