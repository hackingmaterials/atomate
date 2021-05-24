"""
Compute the quasi harmonic approximation gibbs free energy using phonopy.
"""

import sys
import json
from pymongo import MongoClient

import numpy as np

from pymatgen.core import Structure

__author__ = "Kiran Mathew"
__email__ = "kmathew@lbl.gov"

# TODO: @matk86 - I am not sure exactly what this is, but the organization is incorrect here.
# Perhaps in vasp package of atomate? Or pymatgen? -computron


# TODO: @matk86 - please remove these get_db() and get_connection() style methods everywhere.
# They are either already there in pymatgen-db or you can create a single function in common
# utils of atomate. But pretty sure the former. -computron

def get_db(db_file):
    """
    connect to the database and return the connection
    """
    with open(db_file) as f:
        creds = json.loads(f.read())
        conn = MongoClient(creds["host"], creds["port"])
        db = conn[creds["database"]]
        if "admin_user" in creds:
            db.authenticate(creds["admin_user"], creds["admin_password"])
        return db


def get_collection(db_file):
    """
    connect to the database and return task collection
    """
    db = get_db(db_file)
    with open(db_file) as f:
        creds = json.loads(f.read())
        return db[creds["collection"]]


def get_phonopy(structure):
    from phonopy import Phonopy
    from phonopy.structure.atoms import Atoms as PhonopyAtoms

    phon_atoms = PhonopyAtoms(symbols=[str(s.specie) for s in structure],
                              scaled_positions=structure.frac_coords)
    phon_atoms.set_cell(structure.lattice.matrix)
    # supercell size
    scell = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    return Phonopy(phon_atoms, scell)


def get_data(db_file, query):
    coll = get_collection(db_file)
    docs = coll.find(query)
    energies = []
    volumes = []
    force_constants = []
    for d in docs:
        s = Structure.from_dict(d["calcs_reversed"][-1]["output"]['structure'])
        energies.append(d["calcs_reversed"][-1]["output"]['energy'])
        volumes.append(s.volume)
        force_constants.append(d["calcs_reversed"][-1]["output"]['force_constants'])
    return energies, volumes, force_constants


def get_gibbs(structure, db_file, eos="vinet", t_step=10, t_min=0, t_max=1000, mesh=(20, 20, 20),
              plot=False):
    # other eos options: birch_murnaghan, murnaghan
    # The physical units of V and T are \AA^3 and K, respectively.
    # The unit of eV for Helmholtz and Gibbs energies,
    # J/K/mol for C_V and entropy, GPa for for bulk modulus and pressure are used.
    from phonopy import PhonopyQHA

    phonon = get_phonopy(structure)
    energies, volumes, force_constants = get_data(db_file, query={
        "task_label": {"$regex": "gibbs*"}, "formula_pretty": structure.composition.reduced_formula})

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

    phonopy_qha = PhonopyQHA(volumes, energies, eos=eos, temperatures=temperatures[0],
                             free_energy=np.array(free_energy).T, cv=np.array(cv).T,
                             entropy=np.array(entropy).T, t_max=np.max(temperatures[0]))

    # gibbs free energy
    max_t_index = phonopy_qha._qha._max_t_index
    G = phonopy_qha.get_gibbs_temperature()[:max_t_index]
    T = phonopy_qha._qha._temperatures[:max_t_index]
    if plot:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            import matplotlib.pyplot as plt
            plt.plot(T, G)
            plt.savefig("Gibbs.pdf")
            plt.show()
            #phonopy_qha.plot_qha(thin_number=10, volume_temp_exp=None).show()
    else:
        return T, G

# TODO: @matk86 please cleanup, e.g. into an actual unit test -computron
if __name__ == "__main__":
    import os
    from pymatgen.util.testing import PymatgenTest

    structure = PymatgenTest.get_structure("Si")
    db_path = os.path.expanduser("~/db.json")
    get_gibbs(structure, db_path, plot=True)
