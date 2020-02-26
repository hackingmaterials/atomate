# coding: utf-8

import numpy as np
from monty.dev import requires

from pymatgen import Structure
from pymatgen.io.phonopy import get_phonopy_structure
from pymatgen.phonon.bandstructure import (
    PhononBandStructure,
    PhononBandStructureSymmLine,
)
from pymatgen.phonon.dos import CompletePhononDos, PhononDos
from pymatgen.symmetry.bandstructure import HighSymmKpath

try:
    import phonopy
except ImportError:
    phonopy = False

__author__ = "Kiran Mathew, Alex Ganose"
__email__ = "kmathew@lbl.gov, aganose@lbl.gov"

MESH_DENSITY = 100.0  # should always be a float

# TODO: @matk86 - unit tests?


def get_phonopy_gibbs(
    energies,
    volumes,
    force_constants,
    structure,
    t_min,
    t_step,
    t_max,
    mesh,
    eos,
    pressure=0,
):
    """
    Compute QHA gibbs free energy using the phonopy interface.

    Args:
        energies (list):
        volumes (list):
        force_constants (list):
        structure (Structure):
        t_min (float): min temperature
        t_step (float): temperature step
        t_max (float): max temperature
        mesh (list/tuple): reciprocal space density
        eos (str): equation of state used for fitting the energies and the
            volumes. options supported by phonopy: vinet, murnaghan,
            birch_murnaghan
        pressure (float): in GPa, optional.

    Returns:
        (numpy.ndarray, numpy.ndarray): Gibbs free energy, Temperature
    """

    # quasi-harmonic approx
    phonopy_qha = get_phonopy_qha(
        energies,
        volumes,
        force_constants,
        structure,
        t_min,
        t_step,
        t_max,
        mesh,
        eos,
        pressure=pressure,
    )

    # gibbs free energy and temperature
    max_t_index = phonopy_qha._qha._max_t_index
    G = phonopy_qha.get_gibbs_temperature()[:max_t_index]
    T = phonopy_qha._qha._temperatures[:max_t_index]
    return G, T


@requires(phonopy, "phonopy is required to calculate the QHA")
def get_phonopy_qha(
    energies,
    volumes,
    force_constants,
    structure,
    t_min,
    t_step,
    t_max,
    mesh,
    eos,
    pressure=0,
):
    """
    Return phonopy QHA interface.

    Args:
        energies (list):
        volumes (list):
        force_constants (list):
        structure (Structure):
        t_min (float): min temperature
        t_step (float): temperature step
        t_max (float): max temperature
        mesh (list/tuple): reciprocal space density
        eos (str): equation of state used for fitting the energies and the
            volumes. options supported by phonopy: vinet, murnaghan,
            birch_murnaghan
        pressure (float): in GPa, optional.

    Returns:
        PhonopyQHA
    """
    from phonopy import Phonopy
    from phonopy.structure.atoms import Atoms as PhonopyAtoms
    from phonopy import PhonopyQHA
    from phonopy.units import EVAngstromToGPa

    phon_atoms = PhonopyAtoms(
        symbols=[str(s.specie) for s in structure],
        scaled_positions=structure.frac_coords,
        cell=structure.lattice.matrix,
    )
    scell = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    phonon = Phonopy(phon_atoms, scell)
    # compute the required phonon thermal properties
    temperatures = []
    free_energy = []
    entropy = []
    cv = []
    for f in force_constants:
        phonon.set_force_constants(-np.array(f))
        phonon.set_mesh(list(mesh))
        phonon.set_thermal_properties(t_step=t_step, t_min=t_min, t_max=t_max)
        t, g, e, c = phonon.get_thermal_properties()
        temperatures.append(t)
        free_energy.append(g)
        entropy.append(e)
        cv.append(c)

    # add pressure contribution
    energies = (
        np.array(energies) + np.array(volumes) * pressure / EVAngstromToGPa
    )

    # quasi-harmonic approx
    return PhonopyQHA(
        volumes,
        energies,
        eos=eos,
        temperatures=temperatures[0],
        free_energy=np.array(free_energy).T,
        cv=np.array(cv).T,
        entropy=np.array(entropy).T,
        t_max=np.max(temperatures[0]),
    )


def get_phonopy_thermal_expansion(
    energies,
    volumes,
    force_constants,
    structure,
    t_min,
    t_step,
    t_max,
    mesh,
    eos,
    pressure=0,
):
    """
    Compute QHA thermal expansion coefficient using the phonopy interface.

    Args:
        energies (list):
        volumes (list):
        force_constants (list):
        structure (Structure):
        t_min (float): min temperature
        t_step (float): temperature step
        t_max (float): max temperature
        mesh (list/tuple): reciprocal space density
        eos (str): equation of state used for fitting the energies and the
            volumes. options supported by phonopy: vinet, murnaghan,
            birch_murnaghan
        pressure (float): in GPa, optional.

    Returns:
        (numpy.ndarray, numpy.ndarray): thermal expansion coefficient,
        Temperature
    """

    # quasi-harmonic approx
    phonopy_qha = get_phonopy_qha(
        energies,
        volumes,
        force_constants,
        structure,
        t_min,
        t_step,
        t_max,
        mesh,
        eos,
        pressure=pressure,
    )

    # thermal expansion coefficient and temperature
    max_t_index = phonopy_qha._qha._max_t_index
    alpha = phonopy_qha.get_thermal_expansion()[:max_t_index]
    T = phonopy_qha._qha._temperatures[:max_t_index]
    return alpha, T


@requires(phonopy, "phonopy is required to calculate phonon density of states")
def get_phonon_dos(
    structure: Structure,
    supercell_matrix: np.ndarray,
    force_constants: np.ndarray,
    mesh_density: float = MESH_DENSITY,
    num_dos_steps: int = 200,
) -> CompletePhononDos:
    """
    Get a projected phonon density of states.

    Args:
        structure: A structure.
        supercell_matrix: The supercell matrix used to generate the force
            constants.
        force_constants: The force constants in phonopy format.
        mesh_density: The density of the q-point mesh. See the docstring
            for the ``mesh`` argument in Phonopy.init_mesh() for more details.
        num_dos_steps: Number of frequency steps in the energy grid.

    Returns:
        The density of states.
    """
    from phonopy import Phonopy

    structure_phonopy = get_phonopy_structure(structure)
    phonon = Phonopy(structure_phonopy, supercell_matrix=supercell_matrix)
    phonon.set_force_constants(force_constants)
    phonon.run_mesh(
        mesh_density,
        is_mesh_symmetry=False,
        with_eigenvectors=True,
        is_gamma_center=True,
    )

    # get min, max, step frequency
    frequencies = phonon.get_mesh_dict()["frequencies"]
    freq_min = frequencies.min()
    freq_max = frequencies.max()
    freq_pitch = (freq_max - freq_min) / num_dos_steps

    phonon.run_projected_dos(
        freq_min=freq_min, freq_max=freq_max, freq_pitch=freq_pitch
    )

    dos_raw = phonon.projected_dos.get_partial_dos()
    pdoss = {s: dos for s, dos in zip(structure, dos_raw[1])}

    total_dos = PhononDos(dos_raw[0], dos_raw[1].sum(axis=0))
    return CompletePhononDos(structure, total_dos, pdoss)


@requires(phonopy, "phonopy is required to calculate phonon band structures")
def get_phonon_band_structure(
    structure: Structure,
    supercell_matrix: np.ndarray,
    force_constants: np.ndarray,
    mesh_density: float = 100.0,
) -> PhononBandStructure:
    """
    Get a uniform phonon band structure.

    Args:
        structure: A structure.
        supercell_matrix: The supercell matrix used to generate the force
            constants.
        force_constants: The force constants in phonopy format.
        mesh_density: The density of the q-point mesh. See the docstring
            for the ``mesh`` argument in Phonopy.init_mesh() for more details.

    Returns:
        The uniform phonon band structure.
    """
    from phonopy import Phonopy

    structure_phonopy = get_phonopy_structure(structure)
    phonon = Phonopy(structure_phonopy, supercell_matrix=supercell_matrix)
    phonon.set_force_constants(force_constants)
    phonon.run_mesh(mesh_density, is_mesh_symmetry=False, is_gamma_center=True)
    mesh = phonon.get_mesh_dict()

    return PhononBandStructure(
        mesh["qpoints"], mesh["frequencies"], structure.lattice
    )


@requires(phonopy, "phonopy is required to calculate phonon band structures")
def get_line_mode_phonon_band_structure(
    structure: Structure,
    supercell_matrix: np.ndarray,
    force_constants: np.ndarray,
    line_density: float = 20.0,
    symprec: float = 0.01,
) -> PhononBandStructureSymmLine:
    """
    Get a projected phonon density of states.

    Args:
        structure: A structure.
        supercell_matrix: The supercell matrix used to generate the force
            constants.
        force_constants: The force constants in phonopy format.
        line_density: The density along the high symmetry path.
        symprec: Symmetry precision for determining the band structure path.

    Returns:
        The line mode band structure.
    """
    from phonopy import Phonopy

    structure_phonopy = get_phonopy_structure(structure)
    phonon = Phonopy(structure_phonopy, supercell_matrix=supercell_matrix)
    phonon.set_force_constants(force_constants)

    kpath = HighSymmKpath(structure, symprec=symprec)

    kpoints, labels = kpath.get_kpoints(
        line_density=line_density, coords_are_cartesian=False
    )

    phonon.run_qpoints(kpoints)
    frequencies = phonon.qpoints.get_frequencies().T

    labels_dict = dict([(a, k) for a, k in zip(labels, kpoints) if a != ""])

    return PhononBandStructureSymmLine(
        kpoints, frequencies, structure.lattice, labels_dict=labels_dict
    )
