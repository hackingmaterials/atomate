# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
Defines standardized Fireworks that can be chained easily to perform various
sequences of FEFF calculations.
"""

from fireworks import Firework

from pymatgen.io.feff.sets import MPEXAFSSet, MPXANESSet, MPELNESSet

from matmethods.feff.firetasks.write_inputs import WriteFeffFromIOSet
from matmethods.feff.firetasks.run_calc import RunFeffDirect
from matmethods.feff.firetasks.parse_outputs import AbsorptionSpectrumToDbTask
from matmethods.feff.firetasks.run_calc_tscc import RunFeffTscc

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'


class EXAFSFW_tscc(Firework):
    def __init__(self, absorbing_atom, structure, edge="K", radius=10.0, name="EXAFS spectroscopy",
                 feff_input_set=None, feff_cmd="feff", override_default_feff_params=None,
                 db_file=None, parents=None, **kwargs):
        """
        Write the input set for FEFF-EXAFS spectroscopy, run feff and insert the absorption
        coefficient to the database('xas' collection).

        Args:
            absorbing_atom (str): absorbing atom symbol
            structure (Structure): input structure
            radius (float): cluster radius in angstroms
            name (str)
            feff_input_set (FeffDictSet)
            feff_cmd (str): path to the feff binary
            override_default_feff_params (dict): override feff tag settings.
            db_file (str): path to the db file.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """
        override_default_feff_params = override_default_feff_params or {}
        feff_input_set = feff_input_set or MPEXAFSSet(absorbing_atom, structure, edge=edge,
                                                      radius=radius, **override_default_feff_params)

        t = []
        t.append(WriteFeffFromIOSet(absorbing_atom=absorbing_atom, structure=structure,
                                    radius=radius, feff_input_set=feff_input_set))
        t.append(RunFeffTscc())
        t.append(AbsorptionSpectrumToDbTask(absorbing_atom=absorbing_atom, structure=structure,
                                            db_file=db_file, spectrum_type="EXAFS", output_file="xmu.dat"))
        super(EXAFSFW_tscc, self).__init__(t, parents=parents, name="{}-{}".
                                         format(structure.composition.reduced_formula, name), **kwargs)