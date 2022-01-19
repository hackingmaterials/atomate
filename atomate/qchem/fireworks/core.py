# Defines standardized Fireworks that can be chained easily to perform various
# sequences of QChem calculations.

from itertools import chain
import os
import copy

import copy

from fireworks import Firework

from atomate.qchem.firetasks.critic2 import ProcessCritic2, RunCritic2
from atomate.qchem.firetasks.fragmenter import FragmentMolecule
from atomate.qchem.firetasks.parse_outputs import QChemToDb
from atomate.qchem.firetasks.run_calc import RunQChemCustodian
from atomate.qchem.firetasks.write_inputs import WriteInputFromIOSet
from atomate.qchem.firetasks.fragmenter import FragmentMolecule
from atomate.qchem.firetasks.critic2 import RunCritic2, ProcessCritic2
from atomate.qchem.firetasks.geo_transformations import PerturbGeometry

__author__ = "Samuel Blau, Evan Spotte-Smith"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "5/23/18"
__credits__ = "Brandon Wood, Shyam Dwaraknath"


class SinglePointFW(Firework):
    def __init__(
        self,
        molecule=None,
        name="single point",
        qchem_cmd=">>qchem_cmd<<",
        multimode=">>multimode<<",
        max_cores=">>max_cores<<",
        qchem_input_params=None,
        db_file=None,
        parents=None,
        **kwargs
    ):
        """

        Args:
            molecule (Molecule): Input molecule.
            name (str): Name for the Firework.
            qchem_cmd (str): Command to run QChem. Supports env_chk.
            multimode (str): Parallelization scheme, either openmp or mpi. Supports env_chk.
            max_cores (int): Maximum number of cores to parallelize over. Supports env_chk.
            qchem_input_params (dict): Specify kwargs for instantiating the input set parameters.
                Basic uses would be to modify the default inputs of the set, such as dft_rung,
                basis_set, pcm_dielectric, scf_algorithm, or max_scf_cycles. See
                pymatgen/io/qchem/sets.py for default values of all input parameters. For
                instance, if a user wanted to use a more advanced DFT functional, include a pcm
                with a dielectric of 30, and use a larger basis, the user would set
                qchem_input_params = {"dft_rung": 5, "pcm_dielectric": 30, "basis_set":
                "6-311++g**"}. However, more advanced customization of the input is also
                possible through the overwrite_inputs key which allows the user to directly
                modify the rem, pcm, smd, and solvent dictionaries that QChemDictSet passes to
                inputs.py to print an actual input file. For instance, if a user wanted to set
                the sym_ignore flag in the rem section of the input file to true, then they
                would set qchem_input_params = {"overwrite_inputs": "rem": {"sym_ignore":
                "true"}}. Of course, overwrite_inputs could be used in conjunction with more
                typical modifications, as seen in the test_double_FF_opt workflow test.
            db_file (str): Path to file specifying db credentials to place output parsing.
            parents ([Firework]): Parents of this particular Firework.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """

        qchem_input_params = qchem_input_params or {}
        input_file = "mol.qin"
        output_file = "mol.qout"
        t = []
        t.append(
            WriteInputFromIOSet(
                molecule=molecule,
                qchem_input_set="SinglePointSet",
                input_file=input_file,
                qchem_input_params=qchem_input_params,
            )
        )
        t.append(
            RunQChemCustodian(
                qchem_cmd=qchem_cmd,
                multimode=multimode,
                input_file=input_file,
                output_file=output_file,
                max_cores=max_cores,
                job_type="normal",
            )
        )
        t.append(
            QChemToDb(
                db_file=db_file,
                input_file=input_file,
                output_file=output_file,
                additional_fields={"task_label": name},
            )
        )
        super().__init__(t, parents=parents, name=name, **kwargs)


class ForceFW(Firework):
    def __init__(
        self,
        molecule=None,
        name="force calculation",
        qchem_cmd=">>qchem_cmd<<",
        multimode=">>multimode<<",
        max_cores=">>max_cores<<",
        qchem_input_params=None,
        db_file=None,
        parents=None,
        **kwargs
    ):
        """
        Converge the electron density and calculate the atomic forces, aka the gradient.
        Args:
            molecule (Molecule): Input molecule.
            name (str): Name for the Firework.
            qchem_cmd (str): Command to run QChem. Supports env_chk.
            multimode (str): Parallelization scheme, either openmp or mpi. Defaults to openmp.
            max_cores (int): Maximum number of cores to parallelize over. Supports env_chk.
            qchem_input_params (dict): Specify kwargs for instantiating the input set parameters.
                Basic uses would be to modify the default inputs of the set, such as dft_rung,
                basis_set, pcm_dielectric, scf_algorithm, or max_scf_cycles. See
                pymatgen/io/qchem/sets.py for default values of all input parameters. For
                instance, if a user wanted to use a more advanced DFT functional, include a pcm
                with a dielectric of 30, and use a larger basis, the user would set
                qchem_input_params = {"dft_rung": 5, "pcm_dielectric": 30, "basis_set":
                "6-311++g**"}. However, more advanced customization of the input is also
                possible through the overwrite_inputs key which allows the user to directly
                modify the rem, pcm, smd, and solvent dictionaries that QChemDictSet passes to
                inputs.py to print an actual input file. For instance, if a user wanted to set
                the sym_ignore flag in the rem section of the input file to true, then they
                would set qchem_input_params = {"overwrite_inputs": "rem": {"sym_ignore":
                "true"}}. Of course, overwrite_inputs could be used in conjunction with more
                typical modifications, as seen in the test_double_FF_opt workflow test.
            db_file (str): Path to file specifying db credentials to place output parsing.
            parents ([Firework]): Parents of this particular Firework.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """

        qchem_input_params = qchem_input_params or {}
        input_file = "mol.qin"
        output_file = "mol.qout"
        t = []
        t.append(
            WriteInputFromIOSet(
                molecule=molecule,
                qchem_input_set="ForceSet",
                input_file=input_file,
                qchem_input_params=qchem_input_params,
            )
        )
        t.append(
            RunQChemCustodian(
                qchem_cmd=qchem_cmd,
                multimode=multimode,
                input_file=input_file,
                output_file=output_file,
                max_cores=max_cores,
                job_type="normal",
            )
        )
        t.append(
            QChemToDb(
                db_file=db_file,
                input_file=input_file,
                output_file=output_file,
                additional_fields={"task_label": name},
            )
        )
        super(ForceFW, self).__init__(t, parents=parents, name=name, **kwargs)


class OptimizeFW(Firework):
    def __init__(
        self,
        molecule=None,
        name="structure optimization",
        qchem_cmd=">>qchem_cmd<<",
        multimode=">>multimode<<",
        max_cores=">>max_cores<<",
        qchem_input_params=None,
        db_file=None,
        parents=None,
        **kwargs
    ):
        """
        Optimize the given structure.

        Args:
            molecule (Molecule): Input molecule.
            name (str): Name for the Firework.
            qchem_cmd (str): Command to run QChem. Supports env_chk.
            multimode (str): Parallelization scheme, either openmp or mpi. Defaults to openmp.
            max_cores (int): Maximum number of cores to parallelize over. Supports env_chk.
            qchem_input_params (dict): Specify kwargs for instantiating the input set parameters.
                Basic uses would be to modify the default inputs of the set, such as dft_rung,
                basis_set, pcm_dielectric, scf_algorithm, or max_scf_cycles. See
                pymatgen/io/qchem/sets.py for default values of all input parameters. For
                instance, if a user wanted to use a more advanced DFT functional, include a pcm
                with a dielectric of 30, and use a larger basis, the user would set
                qchem_input_params = {"dft_rung": 5, "pcm_dielectric": 30, "basis_set":
                "6-311++g**"}. However, more advanced customization of the input is also
                possible through the overwrite_inputs key which allows the user to directly
                modify the rem, pcm, smd, and solvent dictionaries that QChemDictSet passes to
                inputs.py to print an actual input file. For instance, if a user wanted to set
                the sym_ignore flag in the rem section of the input file to true, then they
                would set qchem_input_params = {"overwrite_inputs": "rem": {"sym_ignore":
                "true"}}. Of course, overwrite_inputs could be used in conjunction with more
                typical modifications, as seen in the test_double_FF_opt workflow test.
            db_file (str): Path to file specifying db credentials to place output parsing.
            parents ([Firework]): Parents of this particular Firework.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """

        qchem_input_params = qchem_input_params or {}
        input_file = "mol.qin"
        output_file = "mol.qout"
        t = []
        t.append(
            WriteInputFromIOSet(
                molecule=molecule,
                qchem_input_set="OptSet",
                input_file=input_file,
                qchem_input_params=qchem_input_params,
            )
        )
        t.append(
            RunQChemCustodian(
                qchem_cmd=qchem_cmd,
                multimode=multimode,
                input_file=input_file,
                output_file=output_file,
                max_cores=max_cores,
                job_type="normal",
            )
        )
        t.append(
            QChemToDb(
                db_file=db_file,
                input_file=input_file,
                output_file=output_file,
                additional_fields={"task_label": name},
            )
        )
        super().__init__(t, parents=parents, name=name, **kwargs)


class TransitionStateFW(Firework):
    def __init__(self,
                 molecule=None,
                 name="transition state structure optimization",
                 qchem_cmd=">>qchem_cmd<<",
                 multimode=">>multimode<<",
                 max_cores=">>max_cores<<",
                 qchem_input_params=None,
                 db_file=None,
                 parents=None,
                 **kwargs):
        """
        Optimize the given molecule to a saddle point of the potential energy surface (transition
        state).
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
        input_file = "mol.qin"
        output_file = "mol.qout"
        t = list()
        t.append(
            WriteInputFromIOSet(
                molecule=molecule,
                qchem_input_set="TransitionStateSet",
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
        super(TransitionStateFW, self).__init__(
            t,
            parents=parents,
            name=name,
            **kwargs)


class FrequencyFW(Firework):
    def __init__(
        self,
        molecule=None,
        name="frequency calculation",
        qchem_cmd=">>qchem_cmd<<",
        multimode=">>multimode<<",
        max_cores=">>max_cores<<",
        qchem_input_params=None,
        db_file=None,
        parents=None,
        **kwargs
    ):
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
        input_file = "mol.qin"
        output_file = "mol.qout"
        t = []
        t.append(
            WriteInputFromIOSet(
                molecule=molecule,
                qchem_input_set="FreqSet",
                input_file=input_file,
                qchem_input_params=qchem_input_params,
            )
        )
        t.append(
            RunQChemCustodian(
                qchem_cmd=qchem_cmd,
                multimode=multimode,
                input_file=input_file,
                output_file=output_file,
                max_cores=max_cores,
                job_type="normal",
            )
        )
        t.append(
            QChemToDb(
                db_file=db_file,
                input_file=input_file,
                output_file=output_file,
                additional_fields={"task_label": name},
            )
        )
        super().__init__(t, parents=parents, name=name, **kwargs)


class PESScanFW(Firework):
    def __init__(self,
                 molecule=None,
                 name="potential energy surface scan",
                 qchem_cmd=">>qchem_cmd<<",
                 multimode=">>multimode<<",
                 max_cores=">>max_cores<<",
                 qchem_input_params=None,
                 scan_variables=None,
                 db_file=None,
                 parents=None,
                 **kwargs):
        """
        Perform a potential energy surface scan by varying bond lengths, angles,
        and/or dihedral angles in a molecule.
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
                                       print an actual input file.
            scan_variables (dict): dict {str: list}, where the key is the type of variable ("stre"
                                   for bond length, "bend" for angle, "tors" for dihedral angle),
                                   and the list contains all of the variable set information
            db_file (str): Path to file specifying db credentials to place output parsing.
            parents ([Firework]): Parents of this particular Firework.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """

        if scan_variables is None:
            raise ValueError("Some variable input must be given! Provide some "
                             "bond, angle, or dihedral angle information.")

        qchem_input_params = qchem_input_params or dict()
        qchem_input_params["scan_variables"] = scan_variables
        input_file = "mol.qin"
        output_file = "mol.qout"
        t = list()

        t.append(
            WriteInputFromIOSet(
                molecule=molecule,
                qchem_input_set="PESScanSet",
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
        super(PESScanFW, self).__init__(
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
                 linked=True,
                 freq_before_opt=False,
                 perturb_geometry=False,
                 mode=None,
                 scale=1.0,
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
                Basic uses would be to modify the default inputs of the set, such as dft_rung,
                basis_set, pcm_dielectric, scf_algorithm, or max_scf_cycles. See
                pymatgen/io/qchem/sets.py for default values of all input parameters. For
                instance, if a user wanted to use a more advanced DFT functional, include a pcm
                with a dielectric of 30, and use a larger basis, the user would set
                qchem_input_params = {"dft_rung": 5, "pcm_dielectric": 30, "basis_set":
                "6-311++g**"}. However, more advanced customization of the input is also
                possible through the overwrite_inputs key which allows the user to directly
                modify the rem, pcm, smd, and solvent dictionaries that QChemDictSet passes to
                inputs.py to print an actual input file. For instance, if a user wanted to set
                the sym_ignore flag in the rem section of the input file to true, then they
                would set qchem_input_params = {"overwrite_inputs": "rem": {"sym_ignore":
                "true"}}. Of course, overwrite_inputs could be used in conjunction with more
                typical modifications, as seen in the test_double_FF_opt workflow test.
            max_iterations (int): Number of perturbation -> optimization -> frequency
                iterations to perform. Defaults to 10.
            max_molecule_perturb_scale (float): The maximum scaled perturbation that can be
                applied to the molecule. Defaults to 0.3.
            freq_before_opt (bool): If True (default False), run a frequency
                calculation before any opt/ts searches to improve understanding
                of the local potential energy surface. Only use this option if
                linked=True.
            perturb_geometry (bool): If True (default False), then modify the input geometry by some
                translation matrix (N x 3, where N is the number of atoms) before optimizing.
            mode (np.ndarray): If not None (default), then perturb the geometry by this matrix.
                This will be ignored if perturb_geometry is False.
            scale (float): Scaling factor for perturbation
            db_file (str): Path to file specifying db credentials to place output parsing.
            parents ([Firework]): Parents of this particular Firework.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """

        qchem_input_params = qchem_input_params or {}
        input_file = "mol.qin"
        output_file = "mol.qout"
        t = []

        if perturb_geometry:
            t.append(PerturbGeometry(
                molecule=molecule,
                mode=mode,
                scale=scale))

            # Make sure that subsequent firetasks use the perturbed Molecule
            molecule = None

        if freq_before_opt:
            t.append(
                WriteInputFromIOSet(
                    molecule=molecule,
                    qchem_input_set="FreqSet",
                    input_file=input_file,
                    qchem_input_params=qchem_input_params))
        else:
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
                linked=linked,
                freq_before_opt=freq_before_opt
            ))
        t.append(
            QChemToDb(
                db_file=db_file,
                input_file=input_file,
                output_file=output_file,
                additional_fields={
                    "task_label": name,
                    "special_run_type": "frequency_flattener",
                    "linked": linked,
                },
            )
        )
        super().__init__(t, parents=parents, name=name, **kwargs)


class FrequencyFlatteningTransitionStateFW(Firework):
    def __init__(self,
                 molecule=None,
                 name="frequency flattening transition state optimization",
                 qchem_cmd=">>qchem_cmd<<",
                 multimode=">>multimode<<",
                 max_cores=">>max_cores<<",
                 qchem_input_params=None,
                 max_iterations=3,
                 max_molecule_perturb_scale=0.3,
                 linked=True,
                 freq_before_opt=True,
                 perturb_geometry=False,
                 mode=None,
                 scale=1,
                 db_file=None,
                 parents=None,
                 **kwargs):
        """
        Iteratively optimize the transition state structure and flatten imaginary frequencies to
        ensure that the resulting structure is a true transition state.
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
                                  iterations to perform. Defaults to 3. Higher numbers are not
                                  recommended, as they rarely lead to improved performance.
            max_molecule_perturb_scale (float): The maximum scaled perturbation that can be
                                                applied to the molecule. Defaults to 0.3.
            linked (bool): If True (default False), the scratch output from one calculation will be passed
                from one calculation to the next, improving convergence behavior.
            freq_before_opt (bool): If True (default False), run a frequency
                calculation before any opt/ts searches to improve understanding
                of the local potential energy surface. Only use this option if
                linked=True.
            perturb_geometry (bool): If True (default False), then modify the input geometry by some
                translation matrix (N x 3, where N is the number of atoms) before optimizing.
            mode (np.ndarray): If not None (default), then perturb the geometry by this matrix.
                This will be ignored if perturb_geometry is False.
            scale (float): Scaling factor for the geometry perturbation.
            db_file (str): Path to file specifying db credentials to place output parsing.
            parents ([Firework]): Parents of this particular Firework.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """

        qchem_input_params = qchem_input_params or {}
        input_file = "mol.qin"
        output_file = "mol.qout"
        runs = list(chain.from_iterable([["ts_" + str(ii), "freq_" + str(ii)]
                                         for ii in range(10)]))
        if freq_before_opt:
            runs.insert(0, "freq_pre")

        t = list()

        if perturb_geometry:
            t.append(PerturbGeometry(
                molecule=molecule,
                mode=mode,
                scale=scale))

            # Make sure that subsequent firetasks use the perturbed Molecule
            molecule = None

        if freq_before_opt:
            t.append(
                WriteInputFromIOSet(
                    molecule=molecule,
                    qchem_input_set="FreqSet",
                    input_file=input_file,
                    qchem_input_params=qchem_input_params))
        else:
            t.append(
                WriteInputFromIOSet(
                    molecule=molecule,
                    qchem_input_set="TransitionStateSet",
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
                transition_state=True,
                linked=linked,
                freq_before_opt=freq_before_opt))
        t.append(
            QChemToDb(
                db_file=db_file,
                input_file=input_file,
                output_file=output_file,
                runs=runs,
                additional_fields={
                    "task_label": name,
                    "special_run_type": "ts_frequency_flattener",
                    "linked": linked
                }))

        super(FrequencyFlatteningTransitionStateFW, self).__init__(
            t,
            parents=parents,
            name=name,
            **kwargs)


class FragmentFW(Firework):
    def __init__(
        self,
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
        **kwargs
    ):
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
                Basic uses would be to modify the default inputs of the set, such as dft_rung,
                basis_set, pcm_dielectric, scf_algorithm, or max_scf_cycles. See
                pymatgen/io/qchem/sets.py for default values of all input parameters. For
                instance, if a user wanted to use a more advanced DFT functional, include a pcm
                with a dielectric of 30, and use a larger basis, the user would set
                qchem_input_params = {"dft_rung": 5, "pcm_dielectric": 30, "basis_set":
                "6-311++g**"}. However, more advanced customization of the input is also
                possible through the overwrite_inputs key which allows the user to directly
                modify the rem, pcm, smd, and solvent dictionaries that QChemDictSet passes to
                inputs.py to print an actual input file. For instance, if a user wanted to set
                the sym_ignore flag in the rem section of the input file to true, then they
                would set qchem_input_params = {"overwrite_inputs": "rem": {"sym_ignore":
                "true"}}. Of course, overwrite_inputs could be used in conjunction with more
                typical modifications, as seen in the test_double_FF_opt workflow test.
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
                check_db=check_db,
            )
        )
        super().__init__(t, parents=parents, name=name, **kwargs)


class CubeAndCritic2FW(Firework):
    def __init__(
        self,
        molecule=None,
        name="cube and critic2",
        qchem_cmd=">>qchem_cmd<<",
        multimode=">>multimode<<",
        max_cores=">>max_cores<<",
        qchem_input_params=None,
        db_file=None,
        parents=None,
        **kwargs
    ):
        """
        Perform a Q-Chem single point calculation in order to generate a cube file of the electron density
        and then analyze the electron density critical points with the Critic2 package.

        Args:
            molecule (Molecule): Input molecule.
            name (str): Name for the Firework.
            qchem_cmd (str): Command to run QChem. Supports env_chk.
            multimode (str): Parallelization scheme, either openmp or mpi. Supports env_chk.
            max_cores (int): Maximum number of cores to parallelize over. Supports env_chk.
            qchem_input_params (dict): Specify kwargs for instantiating the input set parameters.
                Basic uses would be to modify the default inputs of the set, such as dft_rung,
                basis_set, pcm_dielectric, scf_algorithm, or max_scf_cycles. See
                pymatgen/io/qchem/sets.py for default values of all input parameters. For
                instance, if a user wanted to use a more advanced DFT functional, include a pcm
                with a dielectric of 30, and use a larger basis, the user would set
                qchem_input_params = {"dft_rung": 5, "pcm_dielectric": 30, "basis_set":
                "6-311++g**"}. However, more advanced customization of the input is also
                possible through the overwrite_inputs key which allows the user to directly
                modify the rem, pcm, smd, and solvent dictionaries that QChemDictSet passes to
                inputs.py to print an actual input file. For instance, if a user wanted to set
                the sym_ignore flag in the rem section of the input file to true, then they
                would set qchem_input_params = {"overwrite_inputs": "rem": {"sym_ignore":
                "true"}}. Of course, overwrite_inputs could be used in conjunction with more
                typical modifications, as seen in the test_double_FF_opt workflow test.
            db_file (str): Path to file specifying db credentials to place output parsing.
            parents ([Firework]): Parents of this particular Firework.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """

        qchem_input_params = copy.deepcopy(qchem_input_params) or {}
        qchem_input_params["plot_cubes"] = True
        input_file = "mol.qin"
        output_file = "mol.qout"
        t = []
        t.append(
            WriteInputFromIOSet(
                molecule=molecule,
                qchem_input_set="SinglePointSet",
                input_file=input_file,
                qchem_input_params=qchem_input_params,
            )
        )
        t.append(
            RunQChemCustodian(
                qchem_cmd=qchem_cmd,
                multimode=multimode,
                input_file=input_file,
                output_file=output_file,
                max_cores=max_cores,
                job_type="normal",
            )
        )
        t.append(RunCritic2(molecule=molecule, cube_file="dens.0.cube.gz"))
        t.append(ProcessCritic2(molecule=molecule))
        t.append(
            QChemToDb(
                db_file=db_file,
                input_file=input_file,
                output_file=output_file,
                additional_fields={"task_label": name},
            )
        )
        super().__init__(t, parents=parents, name=name, **kwargs)
