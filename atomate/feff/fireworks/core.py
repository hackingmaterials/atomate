# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
Defines standardized Fireworks that can be chained easily to perform various
sequences of FEFF calculations.
"""

from fireworks import Firework

from atomate.utils.utils import load_class
from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.feff.firetasks.glue_tasks import CopyFeffOutputs
from atomate.feff.firetasks.write_inputs import WriteFeffFromIOSet, WriteEXAFSPaths
from atomate.feff.firetasks.run_calc import RunFeffDirect
from atomate.feff.firetasks.parse_outputs import SpectrumToDbTask, AddPathsToFilepadTask

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
             PassCalcLocs(name=name),
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
             PassCalcLocs(name=name),
             SpectrumToDbTask(absorbing_atom=absorbing_atom, structure=structure,
                              db_file=db_file, spectrum_type=spectrum_type, edge=edge,
                              output_file="eels.dat", metadata=metadata)]

        super(EELSFW, self).__init__(t, parents=parents, name="{}-{}".
                                     format(structure.composition.reduced_formula, name), **kwargs)


class EXAFSPathsFW(Firework):
    def __init__(self, absorbing_atom, structure, paths, degeneracies=None, edge="K", radius=10.0,
                 name="EXAFS Paths", feff_input_set=None, feff_cmd="feff",
                 override_default_feff_params=None, parents=None, filepad_file=None, labels=None,
                 metadata=None, **kwargs):
        """
        Write the input set for FEFF-EXAFS spectroscopy with customized scattering paths, run feff,
        and insert the scattering amplitude output files(feffNNNN.dat files) to filepad.

        Args:
            absorbing_atom (str): absorbing atom symbol
            structure (Structure): input structure
            paths (list): list of paths. A path = list of site indices that defines the path legs.
            degeneracies (list): degeneracy of each path.
            edge (str): absorption edge
            radius (float): cluster radius in angstroms
            name (str)
            feff_input_set (FeffDictSet)
            feff_cmd (str): path to the feff binary
            override_default_feff_params (dict): override feff tag settings.
            parents (Firework): Parents of this particular Firework. FW or list of FWS.
            filepad_file (str): path to the filepad config file.
            labels (list): list of label used to tag the files inserted into filepad.
            metadata (dict): meta data
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """
        override_default_feff_params = override_default_feff_params or {}
        override_default_feff_params.update({"user_tag_settings": {"CONTROL": "0 0 0 0 1 1",
                                                                   "PRINT": "0 0 0 1 0 3"}})

        if not feff_input_set:
            fis_cls = load_class("pymatgen.io.feff.sets", "MP{}Set".format("EXAFS"))
            feff_input_set = fis_cls(absorbing_atom, structure, edge=edge, radius=radius,
                                     **override_default_feff_params)

        t = [CopyFeffOutputs(calc_loc=True),
             WriteFeffFromIOSet(absorbing_atom=absorbing_atom, structure=structure, radius=radius,
                                feff_input_set=feff_input_set),
             WriteEXAFSPaths(feff_input_set=feff_input_set, paths=paths, degeneracies=degeneracies),
             RunFeffDirect(feff_cmd=feff_cmd),
             AddPathsToFilepadTask(filepad_file=filepad_file, labels=labels, metadata=metadata)]

        super(EXAFSPathsFW, self).__init__(t, parents=parents, name="{}-{}".format(
            structure.composition.reduced_formula, name), **kwargs)
