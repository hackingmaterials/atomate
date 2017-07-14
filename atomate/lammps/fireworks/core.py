# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
Defines fireworks to be incorporated into workflows
"""

import os

from fireworks import Firework
from atomate.lammps.firetasks.run_calc import RunLammpsDirect
from atomate.common.firetasks.parse_outputs import ToDbTask
from atomate.lammps.firetasks.write_inputs import WriteLammpsFromIOSet
from atomate.lammps.drones import LammpsForceFieldDrone

__author__ = "Brandon Wood"
__email__ = "b.wood@berkeley.edu"

cwd = os.getcwd()

class BasicFW(Firework):
    def __init__(self, job_name, lammps_input, lammps_data, lammps_cmd,
                 data_filename="lammps.data", user_lammps_settings={}, is_forcefield=False,
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

        t = list()
        t.append(WriteLammpsFromIOSet(job_name=job_name, lammps_input=lammps_input,
                                      lammps_data=lammps_data, is_forcefield=is_forcefield))
        t.append(RunLammpsDirect(lammps_cmd=lammps_cmd))
        t.append(ToDbTask(drone=LammpsForceFieldDrone(), mmdb="atomate.utils.database.CalcDb",
                          db_file=db_file, additional_fields={"task_label": job_name}))
        super(BasicFW, self).__init__(t, parents=parents, name=job_name, **kwargs)
