
from six import string_types

def get_all_absorbing_atoms(absorbing_atom, structure):
    if isinstance(absorbing_atom, string_types):
        ab_atoms = structure.indices_from_symbol(absorbing_atom)
    else:
        ab_atoms = [absorbing_atom]
    return ab_atoms
