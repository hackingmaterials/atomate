# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

from typing import Optional, Tuple, List

import numpy as np

from pymatgen import Structure
from pymatgen.transformations.standard_transformations import PerturbStructureTransformation

__author__ = 'Rees Chang, Alex Ganose'
__email__ = 'rc564@cornell.edu, aganose@lbl.gov'


def generate_perturbed_supercells(
    supercell: Structure,
    min_displacement: float = 0.01,
    max_displacement: float = 0.1,
    num_displacements: int = 10,
    n_cells_per_displacement_distance: int = 1,
    min_random_distance: Optional[float] = None
) -> Tuple[List[Structure], np.ndarray]:
    """
    Generate a list of supercells with perturbed atomic sites.

    Args:
        supercell (Structure): Original structure whose atomic sites are to be
            perturbed.
        max_displacement: Maximum displacement distance for perturbing the
            sites in Angstroms.
        min_displacement: Minimum displacement distance for perturbing the
            sites in Angstroms.
        num_displacements: Number of unique displacement distances to
            try, uniformly distributed between ``min_displacement`` and
            ``max_displacement``.
        n_cells_per_displacement_distance: Number of perturbed
            structures to generate for each unique displacement distance.
        min_random_distance: If None (default), then for a given perturbed
            structure, all atoms will move the same distance from their original
            locations. If float, then for a given perturbed structure, the
            distances the atoms move will be uniformly distributed from a
            minimum distance of ``min_random_distance`` to one of the
            displacement distances uniformly sampled between
            ``min_displacement`` and ``max_displacement``.

    Returns:
        A tuple of the ``(structures, displacements)``. Where ``structures`` is
        a list of randomly displaced structures, and ``displacements`` is
        the displacement distance for each structure.
    """
    dists = np.linspace(min_displacement, max_displacement, num_displacements)
    displacements = np.repeat(dists, n_cells_per_displacement_distance)

    perturbed_supercells = []
    for displacement in displacements:
        pst = PerturbStructureTransformation(displacement, min_random_distance)
        perturbed_supercells.append(pst.apply_transformation(supercell))

    return perturbed_supercells, displacements


def csld_main(options, settings):
    """
    Runs CSLD minimization.

    Changes from original version:
        - Made 'prim' an argument in 'phonon_step()'
        - Moved execution files to this main() function to be called from
          atomate
        - Rewrote 'add_common_parameter' in 'common_main' to treat 'options' as

          a dictionary instead of ArgumentParser
    """
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib.backends.backend_pdf import PdfPages
    import atexit
    import csld
    from csld.symmetry_structure import SymmetrizedStructure

    from csld.lattice_dynamics import init_ld_model

    from csld.common_main import upon_exit, \
        init_training
    from csld.csld_main_functions import phonon_step, \
        save_pot, predict, fit_data


    freq_matrix = None #Rees
    pdfout = PdfPages(
        options['pdfout'].strip()) if options['pdfout'].strip() else None
    atexit.register(upon_exit, pdfout)

    prim = SymmetrizedStructure.init_structure(settings['structure'],
                                               options['symm_step'],
                                               options['symm_prim'],
                                               options['log_level'])
    model = init_ld_model(prim, settings['model'], settings[
        'LDFF'] if 'LDFF' in settings.sections() else {}, options['clus_step'],
                          options['symC_step'], options['ldff_step'])
    Amat, fval = init_training(model, settings['training'], options['train_step'])
    ibest, solutions, rel_err = fit_data(model, Amat, fval, settings['fitting'],
                                options['fit_step'], pdfout)
    if settings.has_section('phonon'):
        phonon, freq_matrix = phonon_step(model, solutions, settings['phonon'],
                             options['phonon_step'], pdfout, prim, return_eigen=True)
    if settings.has_section('export_potential'):
        save_pot(model, solutions[ibest], settings['export_potential'],
                 options['save_pot_step'], phonon)
    if settings.has_section('prediction'):
        predict(model, solutions, settings['prediction'], options['pred_step'])

    #OUTPUT
    # freq_matrix is (nbands, nkpoints) = frequencies. Check for negative entries
    # rel_err is cross validation error in percent
    return rel_err, freq_matrix  # also want to return force constants


default_csld_settings = {
    "structure": {"sym_tol": "1e-3"},
    "model": {
        'model_type': 'LD',
        'cluster_in': 'clusters.out',
        'cluster_out': 'clusters.out',
        'symC_in': 'Cmat.mtx',
        'symC_out': 'Cmat.mtx',
        'fractional_distance': False,
    },
    "training":{
        'interface': 'VASP',
        'corr_type': 'f',
        'corr_in': 'Amat.mtx',
        'corr_out': 'Amat.mtx',
        'fval_in': 'fval.txt',
        'fval_out': 'fval.txt'
    },
    "fitting": {
        'solution_in': 'solution_all',
        'solution_out': 'solution_all',
        'nsubset': 5,
        'holdsize': 0.1,
        ## 1 FPC 2 FPC sparse 3 split 4 sparse split
        ## 5 split+ right preconditioning 6 sparse split + r preconditioning
        ## 101 Bayesian CS
        'method': 5,
        # For weight of L1 or L2 regularization
        'mulist': '1E-5 1E-6 1E-7 1E-9',
        'maxIter': 300,
        'tolerance': 1E-6,
        'subsetsize': 0.85,
        'lambda': 0.5,
        'uscale_list': '0.03',
    },
    "phonon": {
        'qpoint_fractional': False,
        # 'Auto' or something like 
        # "[[10,  [0,0,0],'\\Gamma', [0.5,0.5,0.5], 'X', [0.5,0.5,0], 'K']]"
        'wavevector': 'Auto',
        'unit': 'cm',  # THz, meV, eV, cm
        'dos_grid': '15 15 15',  # number of grid points
        'nE_dos': 500, # number of points in DOS
        'ismear': -1,  # 0 (Gaussian), 1 (Lorentzian), -1 (tetrahedron method)
        'epsilon': 0.05,  # width in THz of Gaussian/Lorentzian smearing
        'pdos': True,
        'thermal_T_range': '50 800 50',
        'thermal_out': 'thermal_out.txt'
    },
   'prediction': {
        'interface': 'VASP',
        'corr_type': 'f',
        'corr_in': 'Amat_pred.mtx',
        'corr_out': 'Amat_pred.mtx',
        'fval_in': 'fval_pred.txt',
        'fval_out': 'fval_pred.txt',
        'traindat0': 'fcc222/POSCAR fcc222/traj*'
    },
    'export_potential': {},
    'DEFAULT': {
        "qpoint_fractional": False,
        "true_v_fit": 'true_fit.txt',
        'epsilon': '0.05',
        "bcs_reweight": 'True',
        "bcs_penalty": 'arctan',
        "bcs_jcutoff": '1E-8'
    }
}

default_csld_options = {
    'pdfout': 'plots.pdf',
    'ldff_step': 0,
    'phonon_step': 1,
    'phonon': False,
    'save_pot_step': 1,
    'pot': False,
    'log_level': 1,
    'symm_step': 2,
    'symm_prim': True,
    'clus_step': 3,
    'symC_step': 3,
    'train_step': 3,
    'fit_step': 3,
    'pred_step': 0,
    'refit': False,
    'cont': False,
    'predict': False
}