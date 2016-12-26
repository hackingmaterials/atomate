# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
Defines standardized Fireworks that can be chained easily to perform various
sequences of FEFF calculations.
"""

from fireworks import Firework

from atomate.utils.utils import load_class
from atomate.feff.firetasks.write_inputs import WriteFeffFromIOSet
from atomate.feff.firetasks.run_calc import RunFeffDirect
from atomate.feff.firetasks.parse_outputs import SpectrumToDbTask

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'


class XASFW(Firework):
    def __init__(self, absorbing_atom, structure, spectrum_type, edge="K", radius=10.0,
                 name="XAS spectroscopy", feff_input_set=None, feff_cmd="feff",
                 override_default_feff_params=None, db_file=None, parents=None, metadata=None,
                 **kwargs):
        """
        Write the input set for FEFF-XAS spectroscopy, run feff and insert the absorption
        coefficient to the database(or dump to a json file if db_file=None).

        Args:
            absorbing_atom (str): absorbing atom symbol
            structure (Structure): input structure
            spectrum_type (str): "EXAFS" or "XANES"
            edge (str): absorption edge
            radius (float): cluster radius in angstroms
            name (str)
            feff_input_set (FeffDictSet)
            feff_cmd (str): path to the feff binary
            override_default_feff_params (dict): override feff tag settings.
            db_file (str): path to the db file.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            metadata (dict): meta data
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """
        override_default_feff_params = override_default_feff_params or {}

        if not feff_input_set:
            fis_cls = load_class("pymatgen.io.feff.sets", "MP{}Set".format(spectrum_type))
            feff_input_set = fis_cls(absorbing_atom, structure, edge=edge, radius=radius,
                                     **override_default_feff_params)

        t = [WriteFeffFromIOSet(absorbing_atom=absorbing_atom, structure=structure, radius=radius,
                                feff_input_set=feff_input_set),

             RunFeffDirect(feff_cmd=feff_cmd),

             SpectrumToDbTask(absorbing_atom=absorbing_atom, structure=structure,
                              db_file=db_file, spectrum_type=spectrum_type, edge=edge,
                              output_file="xmu.dat", metadata=metadata)]

        super(XASFW, self).__init__(t, parents=parents, name="{}-{}".
                                    format(structure.composition.reduced_formula, name), **kwargs)


class EELSFW(Firework):
    def __init__(self, absorbing_atom, structure, spectrum_type, edge="K", radius=10.,
                 name="EELS spectroscopy", beam_energy=100, beam_direction=None, collection_angle=1,
                 convergence_angle=1, user_eels_settings=None, feff_input_set=None, feff_cmd="feff",
                 override_default_feff_params=None, db_file=None, parents=None, metadata=None,
                 **kwargs):
        """
        Write the input set for FEFF-EELSS spectroscopy, run feff and insert the core-loss spectrum
        to the database(or dump to a json file if db_file=None).

        Args:
            absorbing_atom (str): absorbing atom symbol
            structure (Structure): input structure
            spectrum_type (str): "ELNES" or "EXELFS"
            edge (str): absorption edge
            radius (float): cluster radius in angstroms
            name (str)
            feff_input_set (FeffDictSet)
            feff_cmd (str): path to the feff binary
            override_default_feff_params (dict): override feff tag settings.
            db_file (str): path to the db file.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            metadata (dict): meta data
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """
        override_default_feff_params = override_default_feff_params or {}

        if not feff_input_set:
            fis_cls = load_class("pymatgen.io.feff.sets", "MP{}Set".format(spectrum_type))
            feff_input_set = fis_cls(absorbing_atom, structure, edge, radius, beam_energy,
                                     beam_direction, collection_angle, convergence_angle,
                                     user_eels_settings=user_eels_settings,
                                     **override_default_feff_params)

        t = [WriteFeffFromIOSet(absorbing_atom=absorbing_atom, structure=structure, radius=radius,
                                feff_input_set=feff_input_set),

             RunFeffDirect(feff_cmd=feff_cmd),

             SpectrumToDbTask(absorbing_atom=absorbing_atom, structure=structure,
                              db_file=db_file, spectrum_type=spectrum_type, edge=edge,
                              output_file="eels.dat", metadata=metadata)]

        super(EELSFW, self).__init__(t, parents=parents, name="{}-{}".
                                     format(structure.composition.reduced_formula, name), **kwargs)
