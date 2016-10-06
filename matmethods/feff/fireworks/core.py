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

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'


class EXAFSFW(Firework):
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
        t.append(RunFeffDirect(feff_cmd=feff_cmd))
        t.append(AbsorptionSpectrumToDbTask(absorbing_atom=absorbing_atom, structure=structure,
                                            db_file=db_file, spectrum_type="EXAFS", output_file="xmu.dat"))
        super(EXAFSFW, self).__init__(t, parents=parents, name="{}-{}".
                                         format(structure.composition.reduced_formula, name), **kwargs)


class XANESFW(Firework):
    def __init__(self, absorbing_atom, structure, edge="K", radius=10.0, name="XANES spectroscopy",
                 feff_input_set=None, feff_cmd="feff", override_default_feff_params=None,
                 db_file=None, parents=None, **kwargs):
        """
        Write the input set for FEFF-XANES spectroscopy and run feff.

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
        feff_input_set = feff_input_set or MPXANESSet(absorbing_atom, structure, edge=edge,
                                                      radius=radius, **override_default_feff_params)

        t = []
        t.append(WriteFeffFromIOSet(absorbing_atom=absorbing_atom, structure=structure,
                                    radius=radius, feff_input_set=feff_input_set))
        t.append(RunFeffDirect(feff_cmd=feff_cmd))
        t.append(AbsorptionSpectrumToDbTask(absorbing_atom=absorbing_atom, structure=structure,
                                            db_file=db_file, spectrum_type="XANES", output_file="xmu.dat"))
        super(XANESFW, self).__init__(t, parents=parents, name="{}-{}".
                                         format(structure.composition.reduced_formula, name), **kwargs)


class ELNESFW(Firework):
    def __init__(self, absorbing_atom, structure, name="ELNES spectroscopy", edge="K", radius=10.,
                 beam_energy=100, beam_direction=None, collection_angle=1, convergence_angle=1,
                 user_eels_settings=None, feff_input_set=None, feff_cmd="feff",
                 override_default_feff_params=None, db_file=None, parents=None, **kwargs):
        """
        Write the input set for FEFF-ELNES spectroscopy and run feff.

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
        mp_id = None
        if kwargs.get("mp_id", None):
            mp_id = kwargs["mp_id"]
            del kwargs["mp_id"]
        override_default_feff_params = override_default_feff_params or {}
        feff_input_set = feff_input_set or MPELNESSet(absorbing_atom, structure, edge, radius,
                                                      beam_energy, beam_direction, collection_angle,
                                                      convergence_angle,
                                                      user_eels_settings=user_eels_settings,
                                                      **override_default_feff_params)

        t = []
        t.append(WriteFeffFromIOSet(absorbing_atom=absorbing_atom, structure=structure,
                                    radius=radius, feff_input_set=feff_input_set))
        t.append(RunFeffDirect(feff_cmd=feff_cmd))
        t.append(AbsorptionSpectrumToDbTask(absorbing_atom=absorbing_atom, structure=structure,
                                            db_file=db_file, spectrum_type="ELNES",
                                            output_file="eels.dat", mp_id=mp_id))
        super(ELNESFW, self).__init__(t, parents=parents, name="{}-{}".
                                         format(structure.composition.reduced_formula, name), **kwargs)
