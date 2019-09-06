# coding: utf-8

# This module defines the torsion potential workflow


from fireworks import Workflow
from atomate.qchem.fireworks.core import OptimizeFW
from atomate.utils.utils import get_logger
from atomate.qchem.firetasks.geo_transformations import RotateTorsion
from atomate.qchem.firetasks.write_inputs import WriteCustomInput

__author__ = "Brandon Wood"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Brandon Wood"
__email__ = "b.wood@berkeley.edu"
__status__ = "Alpha"
__date__ = "5/20/18"
__credits__ = "Sam Blau, Shyam Dwaraknath"

logger = get_logger(__name__)


def get_wf_torsion_potential(molecule,
                             atom_indexes,
                             angles,
                             rem,
                             name="torsion_potential",
                             qchem_cmd=">>qchem_cmd<<",
                             multimode=">>multimode<<",
                             max_cores=">>max_cores<<",
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
            rem (list of two rem dictionaries): a list with two rem dictionaries, one for the first optimization and
            one for the second constrained optimization
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
            **kwargs: Other kwargs that are passed to Firework.__init__.

    Returns: Workflow
    """
    fws = []

    # Optimize the starting molecule fw1
    fw1 = OptimizeFW(
        molecule=molecule,
        name="initial_opt",
        qchem_cmd=qchem_cmd,
        multimode=multimode,
        max_cores=max_cores,
        db_file=db_file,
        **kwargs)
    for idx_t, t in enumerate(fw1.tasks):
        if "WriteInputFromIOSet" in str(t):
            fw1.tasks[idx_t] = WriteCustomInput(molecule=molecule, rem=rem[0])
    fws.append(fw1)

    # Loop to generate all the different rotated molecule optimizations
    for angle in angles:
        rot_opt_fw = OptimizeFW(
            name=("opt_" + str(int(angle))),
            qchem_cmd=qchem_cmd,
            multimode=multimode,
            max_cores=max_cores,
            db_file=db_file,
            parents=fw1,
            **kwargs)
        rot_task = RotateTorsion(atom_indexes=atom_indexes, angle=angle)
        rot_opt_fw.tasks.insert(0, rot_task)
        # define opt section
        opt_line = "tors {a} {b} {c} {d} {ang}".format(
            a=atom_indexes[0],
            b=atom_indexes[1],
            c=atom_indexes[2],
            d=atom_indexes[3],
            ang=angle)
        opt = {"CONSTRAINT": [opt_line]}
        for idx_t, t in enumerate(rot_opt_fw.tasks):
            if "WriteInputFromIOSet" in str(t):
                rot_opt_fw.tasks[idx_t] = WriteCustomInput(rem=rem[1], opt=opt)
        fws.append(rot_opt_fw)

    wfname = "{}:{}".format(molecule.composition.reduced_formula, name)

    return Workflow(fws, name=wfname)
