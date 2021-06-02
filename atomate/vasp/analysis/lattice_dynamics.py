from itertools import product
from typing import Dict, List, Tuple

import numpy as np
import scipy as sp

from atomate.utils.utils import get_logger
<<<<<<< HEAD
<<<<<<< HEAD
from pymatgen import Structure
=======
=======

from hiphive import (ForceConstants, ForceConstantPotential,
                     enforce_rotational_sum_rules, ClusterSpace,
                     StructureContainer)
>>>>>>> 1e75a82a (imports)
from hiphive.cutoffs import is_cutoff_allowed, estimate_maximum_cutoff
from hiphive.fitting import Optimizer

from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor
>>>>>>> 8f4f7e37 (gruneisen and CTE)
from pymatgen.io.phonopy import get_phonopy_structure

from phonopy import Phonopy
from phono3py.phonon3.gruneisen import Gruneisen


__author__ = "Alex Ganose, Rees Chang, Junsoo Park"
__email__ = "aganose@lbl.gov, rc564@cornell.edu, jsyony37@lbl.gov"

logger = get_logger(__name__)

MAX_IMAGINARY_FREQ = 10  # in THz
IMAGINARY_TOL = 0.025  # in THz
#MAX_N_IMAGINARY = np.inf

T_QHA = [i*100 for i in range(16)]
T_RENORM = [i*100 for i in range(0,16)]
T_KLAT = [i*100 for i in range(0,16)]

FIT_METHOD = "rfe" #"least-squares"
eV2J = 1.602e-19

def get_cutoffs(structure: Structure):
    """
    Get a list of trial cutoffs based on a structure.

    An initial guess for the lower bound of the cutoffs is made based on the
    period of the lightest element in the structure, according to:

    ====== === === ===
    .        Cutoff
    ------ -----------
    Period 2ND 3RD 4TH
    ====== === === ===
     1     5.0 3.0 2.0
     2     5.5 3.5 2.5
     3     6.0 4.0 3.0
     4     6.5 4.5 3.0
     5     7.0 5.0 3.5
     6     8.0 5.5 4.0
     7     9.0 6.0 4.5
    ====== === === ===

    The maximum cutoff for each order is determined by the minimum cutoff and
    the following table. A full grid of all possible cutoff combinations is
    generated based on the step size in the table below

    ====== ==== =====
    Cutoff Max  Step
    ====== ==== =====
    2ND    +3.0 0.5
    3RD    +1.5 0.25
    4TH    +1.0 0.25
    ====== ==== =====

    Args:
        structure: A structure.

    Returns:
        A list of trial cutoffs.
    """
    # indexed as min_cutoffs[order][period]
    min_cutoffs = {
        2: {1: 5.0, 2: 5.5, 3: 6.0, 4: 6.5, 5: 7.0, 6: 8.0, 7: 9.0},
        3: {1: 3.0, 2: 3.5, 3: 4.0, 4: 4.5, 5: 5.0, 6: 5.5, 7: 6.0},
        4: {1: 2.0, 2: 2.5, 3: 3.0, 4: 3.0, 5: 3.5, 6: 4.0, 7: 4.5},
    }
    inc = {2: 3, 3: 1.5, 4: 1}
    steps = {2: 0.5, 3: 0.25, 4: 0.25}

    row = min([s.row for s in structure.species])
    mins = {
        2: min_cutoffs[2][row], 3: min_cutoffs[3][row], 4: min_cutoffs[4][row]
    }

    range_two = np.arange(mins[2], mins[2] + inc[2] + steps[2], steps[2])
    range_three = np.arange(mins[3], mins[3] + inc[3] + steps[3], steps[3])
    range_four = np.arange(mins[4], mins[4] + inc[4] + steps[4], steps[4])

<<<<<<< HEAD
    return list(map(list, product(range_two, range_three, range_four)))
=======
    cutoffs = np.array(list(map(list, product(range_two, range_three, range_four))))
    max_cutoff = estimate_maximum_cutoff(AseAtomsAdaptor.get_atoms(supercell_structure))
    logger.info('MAX_CUTOFF \n {}'.format(max_cutoff))    
    good_cutoffs = np.all(cutoffs < np.around(max_cutoff, 4) - 0.0001, axis=1)
    return cutoffs[good_cutoffs].tolist()
>>>>>>> c094a175 (cutoff vs cell_size)


def fit_force_constants(
    parent_structure: Structure,
    supercell_matrix: np.ndarray,
    structures: List["Atoms"],
    all_cutoffs: List[List[float]],
    bulk_modulus: float,
    imaginary_tol: float = IMAGINARY_TOL,
#    max_n_imaginary: int = MAX_N_IMAGINARY,
    max_imaginary_freq: float = MAX_IMAGINARY_FREQ,
    fit_method: str = FIT_METHOD,
<<<<<<< HEAD
) -> Tuple["SortedForceConstants", Dict]:
=======
    n_jobs: int = -1,
    fit_kwargs: Optional[Dict] = None
) -> Tuple["SortedForceConstants", np.ndarray, ClusterSpace, Dict]:
>>>>>>> fddd9681 (get_cutoffs)
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

    Returns:
        A tuple of the best fitted force constants as a hiphive
        ``SortedForceConstants`` object, array of parameters, cluster space,
        and a dictionary of information on the fitting results.
    """
<<<<<<< HEAD
    from hiphive.fitting import Optimizer
    from hiphive import ForceConstantPotential, enforce_rotational_sum_rules

    logger.info("Starting fitting force constants.")
=======
    logger.info("Starting force constant fitting.")
>>>>>>> 8f4f7e37 (gruneisen and CTE)

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
#        "max_n_imaginary": max_n_imaginary,
        "max_imaginary_freq": max_imaginary_freq,
    }

    supercell_atoms = structures[0]

    best_fit = {
        "n_imaginary": np.inf,
<<<<<<< HEAD
        "free_energy": np.inf,
=======
        "rmse_test": np.inf,
        "cluster_space": None,
>>>>>>> cd00c732 (cluster_space to best_fit)
        "force_constants": None,
        "parameters":None,
        "cutoffs": None,
    }
    n_cutoffs = len(all_cutoffs)
<<<<<<< HEAD
    for i, cutoffs in enumerate(all_cutoffs):
=======

    fit_kwargs = fit_kwargs if fit_kwargs else {}
    if fit_method == "rfe" and n_jobs == -1:
        fit_kwargs["n_jobs"] = 1

    cutoff_results = Parallel(n_jobs=5, backend="multiprocessing")(delayed(_run_cutoffs)(
        i, cutoffs, n_cutoffs, parent_structure, structures, supercell_matrix, bulk_modulus,
        fit_method, imaginary_tol, fit_kwargs) for i, cutoffs in enumerate(all_cutoffs))

    logger.info('CUTOFF RESULTS \n {}'.format(cutoff_results))
    
    for result in cutoff_results:
        if result is None:
            continue

        fitting_data["cutoffs"].append(result["cutoffs"])
        fitting_data["rmse_test"].append(result["rmse_test"])
        fitting_data["n_imaginary"].append(result["n_imaginary"])
        fitting_data["min_frequency"].append(result["min_frequency"])
        fitting_data["temperature"].append(result["temperature"])
        fitting_data["free_energy"].append(result["free_energy"])
        fitting_data["entropy"].append(result["entropy"])
        fitting_data["heat_capcity"].append(result["heat_capacity"])

        if (
            result["rmse_test"] < best_fit["rmse_test"]
#            and result["min_frequency"] > -np.abs(max_imaginary_freq)
#            and result["n_imaginary"] <= max_n_imaginary
#            and result["n_imaginary"] < best_fit["n_imaginary"]
        ):
            best_fit.update(result)
            fitting_data["best"] = result["cutoffs"]

    logger.info("Finished fitting force constants.")

    return best_fit["force_constants"], best_fit["parameters"], best_fit["cluster_space"], fitting_data


def _run_cutoffs(
    i,
    cutoffs,
    n_cutoffs,
    parent_structure,
    structures,
    supercell_matrix,
    bulk_modulus,
    fit_method,
    imaginary_tol,
    fit_kwargs
):

    logger.info(
        "Testing cutoffs {} out of {}: {}".format(i + 1, n_cutoffs, cutoffs)
    )
    supercell_atoms = structures[0]

    if not is_cutoff_allowed(supercell_atoms, max(cutoffs)):
>>>>>>> 8f4f7e37 (gruneisen and CTE)
        logger.info(
            "Testing cutoffs {} out of {}: {}".format(i, n_cutoffs, cutoffs)
        )
<<<<<<< HEAD
        sc = get_structure_container(cutoffs, structures)
<<<<<<< HEAD
        opt = Optimizer(sc.get_fit_data(), fit_method)
=======
        logger.info('{}, {}: FINE UNTIL SC'.format(i,cutoffs))
        opt = Optimizer(sc.get_fit_data(), fit_method, **fit_kwargs)
        logger.info('{}, {}: FINE UNTIL OPT'.format(i,cutoffs))
>>>>>>> fddd9681 (get_cutoffs)
=======
        return None

    try:
        sc, cs = get_structure_container(cutoffs, structures)
        ncut = cs.n_dofs
        logger.info('SC and CS generated for cutoff: {}, {}'.format(i,cutoffs))
        opt = Optimizer(sc.get_fit_data(),
                        fit_method,
                        [0,ncut],
                        **fit_kwargs)
        logger.info('Optimizer set up for cutoff: {}, {}'.format(i,cutoffs))
>>>>>>> 79348b89 (bulk_mod input)
        opt.train()
        logger.info('Training complete for cutoff: {}, {}'.format(i,cutoffs))

        parameters = enforce_rotational_sum_rules(
            sc.cluster_space, opt.parameters, ["Huang", "Born-Huang"]
        )
        fcp = ForceConstantPotential(sc.cluster_space, parameters)
        fcs = fcp.get_force_constants(supercell_atoms)
        logger.info('FCS generated for cutoff {}, {}'.format(i,cutoffs))
        n_imaginary, min_freq, free_energy, entropy, Cv, grun, cte, dLfrac = evaluate_force_constants(
            parent_structure, supercell_matrix, fcs, bulk_modulus, T_QHA, imaginary_tol
        )
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD

        fitting_data["cutoffs"].append(cutoffs)
        fitting_data["rmse_test"].append(opt.rmse_test)
        fitting_data["n_imaginary"].append(n_imaginary)
        fitting_data["min_frequency"].append(min_freq)
        fitting_data["300K_free_energy"].append(free_energy)

        if (
            min_freq > -np.abs(max_imaginary_freq)
            and n_imaginary <= max_n_imaginary
            and n_imaginary < best_fit["n_imaginary"]
            and free_energy < best_fit["free_energy"]
        ):
            best_fit.update(
                {
                    "n_imaginary": n_imaginary,
                    "free_energy": free_energy,
                    "force_constants": fcs,
                    "cutoffs": cutoffs,
                }
            )
            fitting_data["best"] = cutoffs

    logger.info("Finished fitting force constants.")

    return best_fit["force_constants"], fitting_data
=======
=======
        logger.info('evaluate_force_constants success!')
>>>>>>> fddd9681 (get_cutoffs)
=======
        logger.info('evaluate_force_constants success: {}, {}'.format(i,cutoffs))
>>>>>>> 79348b89 (bulk_mod input)
        return {
            "cutoffs": cutoffs,
            "rmse_test": opt.rmse_test,
            "n_imaginary": n_imaginary,
            "min_frequency": min_freq,
            "temperature": T_QHA,
            "free_energy": free_energy,
            "entropy": entropy,
            "heat_capacity": Cv,
            "thermal_expansion": cte,
            "cluster_space": sc.cluster_space,
            "parameters": parameters,
            "force_constants": fcs
        }
    except Exception:
        return None
>>>>>>> 8f4f7e37 (gruneisen and CTE)


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

    cs = ClusterSpace(structures[0], cutoffs, symprec=1e-3, acoustic_sum_rules=True)
    logger.debug(cs.__repr__())

    sc = StructureContainer(cs)
    for structure in structures:
        sc.add_structure(structure)
    logger.debug(sc.__repr__())

    return sc, cs


def evaluate_force_constants(
    structure: Structure,
    supercell_matrix: np.ndarray,
    fcs: ForceConstants,
    bulk_modulus: float,
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

    logger.info('EVALUATE FC STARTING')
    fcs2 = fcs.get_fc_array(2)
    fcs3 = fcs.get_fc_array(3)
    parent_phonopy = get_phonopy_structure(structure)
    phonopy = Phonopy(parent_phonopy, supercell_matrix=supercell_matrix)
    vol = phonopy.primitive.get_volume()
    natom = phonopy.primitive.get_number_of_atoms()
    mesh = supercell_matrix.diagonal()*5
    
    phonopy.set_force_constants(fcs2)
    logger.info('FCS2 SET')
    phonopy.set_mesh(mesh,is_eigenvectors=True,is_mesh_symmetry=False) #run_mesh(is_gamma_center=True)
    logger.info('MESH SET')
    phonopy.run_thermal_properties(temperatures=T)
    logger.info('Thermal properties successfully run!')

    _, free_energy, entropy, Cv = phonopy.get_thermal_properties()
    free_energy *= 1000/sp.constants.Avogadro/eV2J/natom # kJ/mol to eV/atom
    entropy *= 1/sp.constants.Avogadro/eV2J/natom # J/K/mol to eV/K/atom
    Cv *= 1/sp.constants.Avogadro/eV2J/natom # J/K/mol to eV/K/atom
    logger.info('Units converted!')
    
    freq = phonopy.mesh.frequencies # in THz
    # find imaginary modes at gamma
#    phonopy.run_qpoints([0, 0, 0])
#    gamma_eigs = phonopy.get_qpoints_dict()["frequencies"]
    logger.info('FREQ \n {}:'.format(freq))
    n_imaginary = int(np.sum(freq < -np.abs(imaginary_tol)))
    logger.info('N_IMAGINARY: {}'.format(n_imaginary))
    min_freq = np.min(freq)

    if n_imeginary > 0 and bulk_modulus is not None:
        logger.info('No imaginary modes! Computing Gruneisen and thermal expansion.')
        grun, cte = gruneisen(phonopy,fcs2,fcs3,mesh,T,Cv,bulk_modulus,vol)
        dLfrac = thermal_expansion(T,cte)
    elif n_imeginary > 0 and bulk_modulus is None: # need bulk modulus input
        logger.warning('No imaginary modes, but bulk modulus is not supplied!')
        logger.warning('Gruneisen and thermal expansion are not calculated.')
        grun = np.zeros((len(T),3))
        cte = np.zeros((len(T),3))
        dLfrac = np.zeros((len(T),3))
    else: # do not calculate these if imaginary modes exist
        logger.warning('Imaginary modes found!')
        grun = np.zeros((len(T),3))
        cte = np.zeros((len(T),3))
        dLfrac = np.zeros((len(T),3))
        
    return n_imaginary, min_freq, free_energy, entropy, Cv, grun, cte, dLfrac


def get_total_grun(
        omega: np.ndarray,
        grun: np.ndarray,
        kweight: np.ndarray,
        T: float
) -> np.ndarray:
    total = 0
    weight = 0
    nptk = omega.shape[0]
    nbands = omega.shape[1]
    omega = abs(omega)*1e12*2*np.pi
    if T==0:
        total = np.zeros((3,3))
        grun_total_diag = np.zeros(3)
    else:
        for i in range(nptk):
            for j in range(nbands):
                x = hbar*omega[i,j]/(2.0*kB*T)
                dBE = (x/np.sinh(x))**2
                weight += dBE*kweight[i]
                total += dBE*kweight[i]*grun[i,j]
        total = total/weight
        grun_total_diag = np.array([total[0,2],total[1,1],total[2,0]])

        def percent_diff(a,b):
            return abs((a-b)/b)
        # This process preserves cell symmetry upon thermal expansion, i.e., it prevents
        # symmetry-identical directions from inadvertently expanding by different ratios
        # when the Gruneisen routine returns slighlty different ratios for those directions
        if percent_diff(grun_total_diag[0],grun_total_diag[1]) < 0.1:
            avg = np.mean((grun_total_diag[0],grun_total_diag[1]))
            grun_total_diag[0] = avg
            grun_total_diag[1] = avg
        elif percent_diff(grun_total_diag[0],grun_total_diag[2]) < 0.1:
            avg = np.mean((grun_total_diag[0],grun_total_diag[2]))
            grun_total_diag[0] = avg
            grun_total_diag[2] = avg
        elif percent_diff(grun_total_diag[1],grun_total_diag[2]) < 0.1:
            avg = np.mean((grun_total_diag[1],grun_total_diag[2]))
            grun_total_diag[1] = avg
            grun_total_diag[2] = avg
        else:
            pass
    return grun_total_diag


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

