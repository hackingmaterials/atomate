# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
Defines fireworks to be incorporated into workflows
"""

from fireworks import Firework
from atomate.lammps.firetasks.run_calc import RunLammpsDirect
from atomate.lammps.firetasks.parse_outputs import LammpsToDBTask
from pymatgen.io.lammps.input import DictLammpsInput

__author__ = "Brandon Wood"
__email__ = "b.wood@berkeley.edu"


class ReadRunFW(Firework):
    def __init__(self, job_name, lammps_input, lammps_data, data_filename,
                 user_lammps_settings={}, is_forcefield=False, lammps_cmd="lammps",
                 db_file=None, parents=None, **kwargs):
        """
        Read, run, and store a lammps simulation. This is useful for lammps simulations that already have
        input and data files written.

        Args:
            job_name: descriptive name for lammps simulation
            lammps_input: path to lammps style input file
            lammps_data: path to lammps data file
            data_filename: data file name
            user_lammps_settings: settings that will overwrite input files
            is_forcefield: whether the data file has forcefield info in it.
                This is required only if lammps_data is a path to the data file instead of a data object
            lammps_cmd: command to run lammps
            db_file: path to file specifying db credentials to place output parsing
            parents ([Fireworks)]: parents of this particular Firework
            \*\*kwargs: other kwargs that are passed to Firework.__init__.
        """

        lammps_cmd = lammps_cmd + " -in " + lammps_input
        lammps_input_set = DictLammpsInput.from_file(job_name, lammps_input, lammps_data,
                                                     data_filename, user_lammps_settings, is_forcefield)
        t = []
        t.append(RunLammpsDirect(lammps_cmd=lammps_cmd))
        t.append(LammpsToDBTask(lammps_input=lammps_input_set, db_file=db_file))
        super(ReadRunFW, self).__init__(t, parents=parents, name=job_name, **kwargs)
