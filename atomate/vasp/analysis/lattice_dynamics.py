from itertools import product
from joblib import Parallel, delayed
from typing import Dict, List, Tuple, Optional

import numpy as np
import scipy as sp

from atomate.utils.utils import get_logger
from hiphive.cutoffs import is_cutoff_allowed, estimate_maximum_cutoff
from hiphive import ForceConstants
from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.phonopy import get_phonopy_structure

__author__ = "Alex Ganose, Rees Chang, Junsoo Park"
__email__ = "aganose@lbl.gov, rc564@cornell.edu, jsyony37@lbl.gov"

logger = get_logger(__name__)

MAX_IMAGINARY_FREQ = 10  # in THz
IMAGINARY_TOL = 0.025  # in THz
MAX_N_IMAGINARY = np.inf
FIT_METHOD = "least-squares"
eV2J = 1.602e-19

def get_cutoffs(supercell_structure: Structure):
    """
    Get a list of trial cutoffs based on a supercell structure.

    An initial guess for the lower bound of the cutoffs is made based on the
    period of the lightest element in the structure, according to:

    ====== === === ===
    .        Cutoff
    ------ -----------
    Period 2ND 3RD 4TH
    ====== === === ===
     1     5.5 3.0 3.0
     2     5.5 3.0 3.0
     3     6.5 4.0 3.5
     4     7.5 5.0 4.0
     5     8.5 6.0 4.5
     6     9.5 7.0 5.0
     7     9.5 7.0 5.0
    ====== === === ===

    The maximum cutoff for each order is determined by the minimum cutoff and
    the following table. A full grid of all possible cutoff combinations is
    generated based on the step size in the table below

    ====== ==== =====
    Cutoff Max  Step
    ====== ==== =====
    2ND    +3.0 0.5
    3RD    +1.5 0.5
    4TH    +1.0 0.5
    ====== ==== =====

    Finally, the max cutoff size is determined by the supercell lattice dimensions.
    Cutoffs which result in multiple of the same orbits being populated will be
    discounted.

    Args:
        supercell_structure: A structure.

    Returns:
        A list of trial cutoffs.
    """
    # indexed as min_cutoffs[order][period]
    min_cutoffs = {
        2: {1: 5.5, 2: 5.5, 3: 6.5, 4: 7.5, 5: 8.5, 6: 9.5, 7: 9.5},
        3: {1: 3.0, 2: 3.0, 3: 4.0, 4: 5.0, 5: 6.0, 6: 7.0, 7: 7.0},
        4: {1: 3.0, 2: 3.0, 3: 3.5, 4: 4.0, 5: 4.5, 6: 5.0, 7: 5.0},
    }
    inc = {2: 3, 3: 1.5, 4: 1}
    steps = {2: 0.5, 3: 0.5, 4: 0.5}

    row = min([s.row for s in supercell_structure.species])
    mins = {
        2: min_cutoffs[2][row], 3: min_cutoffs[3][row], 4: min_cutoffs[4][row]
    }

    range_two = np.arange(mins[2], mins[2] + inc[2] + steps[2], steps[2])
    range_three = np.arange(mins[3], mins[3] + inc[3] + steps[3], steps[3])
    range_four = np.arange(mins[4], mins[4] + inc[4] + steps[4], steps[4])

    cutoffs = np.array(list(map(list, product(range_two, range_three, range_four))))
    max_cutoff = estimate_maximum_cutoff(AseAtomsAdaptor.get_atoms(supercell_structure))
    good_cutoffs = np.all(cutoffs < np.around(max_cutoff, 4) - 0.0001, axis=1)
    return cutoffs[good_cutoffs].tolist()


def fit_force_constants(
    parent_structure: Structure,
    supercell_matrix: np.ndarray,
    structures: List["Atoms"],
    all_cutoffs: List[List[float]],
    imaginary_tol: float = IMAGINARY_TOL,
    max_n_imaginary: int = MAX_N_IMAGINARY,
    max_imaginary_freq: float = MAX_IMAGINARY_FREQ,
    fit_method: str = FIT_METHOD,
    n_jobs: int = -1,
    fit_kwargs: Optional[Dict] = None
) -> Tuple["SortedForceConstants", Dict]:
    """
    Fit force constants using hiphive and select the optimum cutoff values.

    The optimum cutoffs will be determined according to:
    1. Number imaginary modes < ``max_n_imaginary``.
    2. Most negative imaginary frequency < ``max_imaginary_freq``.
    3. Least number of imaginary modes.
    4. Lowest free energy at 300 K.

    If criteria 1 and 2 are not satisfied, None will be returned as the
    force constants.

    Args:
        parent_structure: Initial input structure.
        supercell_matrix: Supercell transformation matrix.
        structures: A list of ase atoms objects with "forces" and
            "displacements" arrays added, as required by hiPhive.
        all_cutoffs: A nested list of cutoff values to trial. Each set of
            cutoffs contains the radii for different orders starting with second
            order.
        imaginary_tol: Tolerance used to decide if a phonon mode is imaginary,
            in THz.
        max_n_imaginary: Maximum number of imaginary modes allowed in the
            the final fitted force constant solution. If this criteria is not
            reached by any cutoff combination this FireTask will fizzle.
        max_imaginary_freq: Maximum allowed imaginary frequency in the
            final fitted force constant solution. If this criteria is not
            reached by any cutoff combination this FireTask will fizzle.
        fit_method: Method used for fitting force constants. This can be
            any of the values allowed by the hiphive ``Optimizer`` class.
        n_jobs: Number of processors to use for fitting coefficients. -1 means use all
            processors.
        fit_kwargs: Additional arguements passed to the hiphive force constant
            optimizer.

    Returns:
        A tuple of the best fitted force constants as a hiphive
        ``SortedForceConstants`` object and a dictionary of information on the
        fitting process.
    """
    logger.info("Starting force constant fitting.")

    fitting_data = {
        "cutoffs": [],
        "rmse_test": [],
        "n_imaginary": [],
        "min_frequency": [],
        "temperature": [],
        "free_energy": [],
        "entropy": [],
        "heat_capacity": [],
        "fit_method": fit_method,
        "imaginary_tol": imaginary_tol,
        "max_n_imaginary": max_n_imaginary,
        "max_imaginary_freq": max_imaginary_freq,
    }

    best_fit = {
        "n_imaginary": np.inf,
        "rmse_test": np.inf,
        "force_constants": None,
        "cutoffs": None,
    }
    n_cutoffs = len(all_cutoffs)

    fit_kwargs = fit_kwargs if fit_kwargs else {}
    if fit_method == "rfe" and n_jobs == -1:
        fit_kwargs["n_jobs"] = 1

    cutoff_results = Parallel(n_jobs=-1, backend="multiprocessing")(delayed(_run_cutoffs)(
        i, cutoffs, n_cutoffs, parent_structure, structures, supercell_matrix,
        fit_method, imaginary_tol, fit_kwargs) for i, cutoffs in enumerate(all_cutoffs))

    for result in cutoff_results:
        if result is None:
            pass

        fitting_data["cutoffs"].append(result["cutoffs"])
        fitting_data["rmse_test"].append(result["rmse_test"])
        fitting_data["n_imaginary"].append(result["n_imaginary"])
        fitting_data["min_frequency"].append(result["min_frequency"])
        fitting_data["temperature"].append(result["temperature"])
        fitting_data["free_energy"].append(result["free_energy"])
        fitting_data["entropy"].append(result["entropy"])
        fitting_data["heat_capcity"].append(result["heat_capacity"])

        if (
            result["min_frequency"] > -np.abs(max_imaginary_freq)
#            and result["n_imaginary"] <= max_n_imaginary
#            and result["n_imaginary"] < best_fit["n_imaginary"]
            and result["rmse_test"] < best_fit["rmse_test"]
        ):
            best_fit.update(result)
            fitting_data["best"] = result["cutoffs"]

    logger.info("Finished fitting force constants.")

    return best_fit["force_constants"], fitting_data


def _run_cutoffs(
    i,
    cutoffs,
    n_cutoffs,
    parent_structure,
    structures,
    supercell_matrix,
    fit_method,
    imaginary_tol,
    fit_kwargs
):
    from hiphive.fitting import Optimizer
    from hiphive import ForceConstantPotential, enforce_rotational_sum_rules

    logger.info(
        "Testing cutoffs {} out of {}: {}".format(i + 1, n_cutoffs, cutoffs)
    )
    supercell_atoms = structures[0]

    if not is_cutoff_allowed(supercell_atoms, max(cutoffs)):
        logger.info(
            "Skipping cutoff due as it is not commensurate with supercell size"
        )
        return None

    try:
        sc = get_structure_container(cutoffs, structures)
        opt = Optimizer(sc.get_fit_data(), fit_method, **fit_kwargs)
        opt.train()

        parameters = enforce_rotational_sum_rules(
            sc.cluster_space, opt.parameters, ["Huang", "Born-Huang"]
        )
        fcp = ForceConstantPotential(sc.cluster_space, parameters)
        fcs = fcp.get_force_constants(supercell_atoms)

        T_qha = [i*100 for i in range(16)]
#        phonopy_fcs = fcs.get_fc_array(order=2)
        n_imaginary, min_freq, free_energy, entropy, Cv, grun, cte, dLfrac = evaluate_force_constants(
            parent_structure, supercell_matrix, fcs, T_qha, imaginary_tol
        )
        return {
            "cutoffs": cutoffs,
            "rmse_test": opt.rmse_test,
            "n_imaginary": n_imaginary,
            "min_frequency": min_freq,
            "temperature": T_qha,
            "free_energy": free_energy,
            "entropy": entropy,
            "heat_capacity": Cv,
            "thermal_expansion": cte,
            "force_constants": fcp
        }
    except Exception:
        return None


def get_structure_container(
        cutoffs: List[float], structures: List["Atoms"]
) -> "StructureContainer":
    """
    Get a hiPhive StructureContainer from cutoffs and a list of atoms objects.

    Args:
        cutoffs: Cutoff radii for different orders starting with second order.
        structures: A list of ase atoms objects with the "forces" and
            "displacements" arrays included.

    Returns:
        A hiPhive StructureContainer.
    """
    from hiphive import ClusterSpace, StructureContainer

    cs = ClusterSpace(structures[0], cutoffs)
    logger.debug(cs.__repr__())

    sc = StructureContainer(cs)
    for structure in structures:
        sc.add_structure(structure)
    logger.debug(sc.__repr__())

    return sc


def evaluate_force_constants(
    structure: Structure,
    supercell_matrix: np.ndarray,
    force_constants: ForceConstants,
    T: List,
    imaginary_tol: float = IMAGINARY_TOL
) -> Tuple[int, float, List, List, List]:
    """
    Uses the force constants to extract phonon properties. Used for comparing
    the accuracy of force constant fits.

    Args:
        structure: The parent structure.
        supercell_matrix: The supercell transformation matrix.
        force_constants: The force constants in numpy format.
        imaginary_tol: Tolerance used to decide if a phonon mode is imaginary,
            in THz.

    Returns:
        A tuple of the number of imaginary modes at Gamma, the minimum phonon
        frequency at Gamma, and the free energy, entropy, and heat capacity
    """
    from phonopy import Phonopy

    fcs2 = fcs.get_fc_array(2)
    fcs3 = fcs.get_fc_array(3)
    parent_phonopy = get_phonopy_structure(structure)
    phonopy = Phonopy(parent_phonopy, supercell_matrix=supercell_matrix)
    vol = phonopy.primitive.get_volume()
    natom = phonopy.primitive.get_number_of_atoms()
    
    phonopy.set_force_constants(fcs2)
    phonopy.run_mesh(is_gamma_center=True)
    phonopy.run_thermal_properties(temperatures=T)

    free_energy = phonopy.get_thermal_properties_dict()["free_energy"]
    free_energy *= 1000/sp.constants.Avogadro/eV2J/natom # kJ/mol to eV/atom
    entropy = phonopy.get_thermal_properties_dict()["entropy"]
    entropy *= 1/sp.constants.Avogadro/eV2J/natom # J/K/mol to eV/K/atom
    Cv = phonopy.get_thermal_properties_dict()["heat_capacity"]
    Cv *= 1/sp.constants.Avogadro/eV2J/natom # J/K/mol to eV/K/atom
    
    freq = phonopy.mesh.frequencies # in THz
    # find imaginary modes at gamma
#    phonopy.run_qpoints([0, 0, 0])
#    gamma_eigs = phonopy.get_qpoints_dict()["frequencies"]
    n_imaginary = int(np.sum(freq < -np.abs(imaginary_tol)))
    min_frequency = np.min(freq)

    if n_imeginary > 0:
        grun, cte = gruneisen(phonopy,fcs2,fcs3,mesh,T,Cv,bulk_mod,vol)
        dLfrac = thermal_expansion(T,cte)
    else: # do not calculate these if imaginary modes exist
        grun = np.zeros((len(T),3))
        cte = np.zeros((len(T),3))
        dLfrac = np.zeros((len(T),3))
        
    return n_imaginary, min_frequency, free_energy, entropy, Cv, grun, cte, dLfrac


def gruneisen(
        phonopy: Phonopy,
        fcs2: np.ndarray,
        fcs3: np.ndarray,
        mesh: List,
        temperature: List,
        Cv: np.ndarray, # in eV/K/atom
        bulk_mod: float, # in GPa
        vol: float # in A^3
) -> Tuple[List,List]:

    from phono3py.phonon3.gruneisen import Gruneisen
    
    gruneisen = Gruneisen(fcs2,fcs3,phonopy.supercell,phonopy.primitive)
    gruneisen.set_sampling_mesh(mesh,is_gamma_center=True)
    gruneisen.run()
    grun = gruneisen.get_gruneisen_parameters() # (nptk,nmode,3,3)
    omega = gruneisen._frequencies
    qp = gruneisen._qpoints
    kweight = gruneisen._weights
    grun_tot = list()
    for temp in temperature:
        grun_tot.append(get_total_grun(omega,grun,kweight,temp))
    grun_tot = np.array(np.nan_to_num(np.array(grun_tot)))
    Cv *= eV2J*phonopy.primitive.get_number_of_atoms() # eV/K/atom to J/K
    # linear thermal expansion coefficient     
    cte = grun_tot*(Cv.repeat(3).reshape((len(Cv),3)))/(vol/10**30)/(bulk_mod*10**9)/3
    cte = np.nan_to_num(cte)    
    logger.info('Frequency : {}'.format(np.sort(omega.flatten())))
    logger.info('Gruneisen : {}'.format(grun_tot))
    logger.info('CTE : {}'.format(cte))    
    return grun_tot, cte


def thermal_expansion(
        temperature: List,
        cte: List,
        T: Optional[float]=None
) -> np.ndarray:
    assert len(temperature)==len(cte)
    if 0 not in temperature:
        temperature = [0] + temperature
        cte = [[0,0,0]] + cte
    temperature = np.array(temperature)
    ind = np.argsort(temperature)
    temperature = temperature[ind]
    cte = np.array(cte)[ind]
    # linear expansion fraction
    dLfrac = copy(cte)
    for t in range(len(temperature)):
        dLfrac[t,:] = np.trapz(cte[:t+1,:],temperature[:t+1],axis=0)
    dLfrac = np.nan_to_num(dLfrac)
    logger.info('dLfrac : {}'.format(dLfrac))
    if T is None:
        return dLfrac
    else:
        try:
            T_ind = np.where(temperature==T)[0][0]
            return np.array([dLfrac[T_ind]])
        except:
            raise ValueError('Designated T does not exist in the temperature array!')

