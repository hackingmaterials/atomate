# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
Defines standardized Fireworks that can be chained easily to perform various
sequences of FEFF calculations.
"""

from fireworks import Firework

from pymatgen.io.feff.sets import MPEXAFSSet, MPXANESSet

from matmethods.feff.firetasks.run_calc import RunFeffDirect
from matmethods.feff.firetasks.write_inputs import WriteFeffFromIOSet

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'


class EXAFSFW(Firework):
    def __init__(self, absorbing_atom, structure, radius=10.0, name="EXAFS spectroscopy",
                 feff_input_set=None, feff_cmd="feff", override_default_feff_params=None,
                 db_file=None, parents=None, **kwargs):
        """

        Args:
            absorbing_atom:
            structure:
            radius:
            name:
            feff_input_set:
            feff_cmd:
            override_default_feff_params:
            db_file:
            parents:
            **kwargs:
        """
        override_default_feff_params = override_default_feff_params or {}
        feff_input_set = feff_input_set or MPEXAFSSet(absorbing_atom, structure, radius,
                                                      **override_default_feff_params)

        t = []
        print(feff_input_set.tags)
        t.append(WriteFeffFromIOSet(absorbing_atom=absorbing_atom, structure=structure,
                                    radius=radius, feff_input_set=feff_input_set))
        t.append(RunFeffDirect(feff_cmd=feff_cmd))
        #t.append(FeffToDbTask(db_file=db_file, additional_fields={"task_label": name}))
        super(EXAFSFW, self).__init__(t, parents=parents, name="{}-{}".
                                         format(structure.composition.reduced_formula, name), **kwargs)


class XANESFW(Firework):
    def __init__(self, absorbing_atom, structure, radius=10.0, name="XANES spectroscopy",
                 feff_input_set=None, feff_cmd="feff", override_default_feff_params=None,
                 db_file=None, parents=None, **kwargs):
        """

        Args:
            absorbing_atom:
            structure:
            radius:
            name:
            feff_input_set:
            feff_cmd:
            override_default_feff_params:
            db_file:
            parents:
            **kwargs:
        """
        override_default_feff_params = override_default_feff_params or {}
        feff_input_set = feff_input_set or MPXANESSet(absorbing_atom, structure, radius,
                                                      **override_default_feff_params)

        t = []
        t.append(WriteFeffFromIOSet(absorbing_atom=absorbing_atom, structure=structure,
                                    radius=radius, feff_input_set=feff_input_set))
        t.append(RunFeffDirect(feff_cmd=feff_cmd))
        #t.append(FeffToDbTask(db_file=db_file, additional_fields={"task_label": name}))
        super(XANESFW, self).__init__(t, parents=parents, name="{}-{}".
                                         format(structure.composition.reduced_formula, name), **kwargs)
