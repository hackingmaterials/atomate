# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
Defines fireworks to be incorporated into workflows
"""

from fireworks import Firework

from atomate.lammps.firetasks.run_calc import RunLammpsDirect
from atomate.common.firetasks.parse_outputs import ToDbTask
from atomate.lammps.firetasks.write_inputs import WriteLammpsFromIOSet
from atomate.lammps.drones import LammpsForceFieldDrone

__author__ = "Brandon Wood"
__email__ = "b.wood@berkeley.edu"


class LammpsFW(Firework):
    def __init__(self, job_name, lammps_input_set, input_filename="lammps.in",
                 data_filename="lammps.data", lammps_cmd="lammps", db_file=None, parents=None, **kwargs):
        """
        Read, run, and store a lammps simulation. This is useful for lammps simulations that already
        have input and data files written.

        Args:
            job_name (str): descriptive name for lammps simulation
            lammps_input (str): path to lammps style input file (either plain or json file).
                The settings erad form the file will be overridden by the user_lammps_settings.
            lammps_data (str/LammpsData): either LammpsData object or path to lammps data file.
            data_filename (str): the name of the file to which the lammps data is to be written.
            user_lammps_settings (dict): settings that will override the dettings read from
                lammps_input.
            is_forcefield (bool): whether the data file has forcefield info in it.
                This is required only if lammps_data is a path to the data file instead of a data
                object.
            lammps_cmd (str): command to run lammps.
            db_file (str): path to file specifying db credentials to place output parsing.
            parents ([Fireworks)]: parents of this particular Firework.
            \*\*kwargs: other kwargs that are passed to Firework.__init__.
        """

        t = [
            WriteLammpsFromIOSet(lammps_input_set=lammps_input_set, input_filename=input_filename,
                                 data_filename=data_filename),
            RunLammpsDirect(lammps_cmd=lammps_cmd),
            ToDbTask(drone=LammpsForceFieldDrone(), mmdb="atomate.lammps.database.LammpsCalcDb",
                     db_file=db_file, additional_fields={"task_label": job_name})
        ]
        super(LammpsFW, self).__init__(t, parents=parents, name=job_name, **kwargs)


# TODO: implement this
class PackmolFW(Firework):
    pass
