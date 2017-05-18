# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

from six import string_types

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'


def get_all_absorbing_atoms(absorbing_atom, structure):
    """
    If the absorbing atom is specified as a string, return the list of all
    sites with the atomic symbol. Else just return back the int index.

    Args:
        absorbing_atom (str/int): site symbol or index number
        structure (Structure)

    Returns:
        list: list of site indices
    """
    # TODO: @matk86 - I'd delete this method and replace with:
    # ab_atoms = [absorbing_atom] if isinstance (absorbing_atom, int) else structure.indices_from_symbol(absorbing_atom)
    # inside the code of wherever it is called. Having a separate method and a separate module
    # (util.py) for just this simple thing makes things complicated. -computron

    if isinstance(absorbing_atom, string_types):
        ab_atoms = structure.indices_from_symbol(absorbing_atom)
    else:
        ab_atoms = [absorbing_atom]
    return ab_atoms
