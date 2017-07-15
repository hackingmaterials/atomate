# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
Defines fireworks to be incorporated into workflows.
"""

from fireworks import Firework

from atomate.lammps.firetasks.run_calc import RunLammpsDirect
from atomate.common.firetasks.parse_outputs import ToDbTask
from atomate.lammps.firetasks.write_inputs import WriteLammpsFromIOSet
from atomate.lammps.drones import LammpsForceFieldDrone

__author__ = "Brandon Wood"
__email__ = "b.wood@berkeley.edu"


class LammpsFW(Firework):
    def __init__(self, lammps_input_set, input_filename="lammps.in", data_filename="lammps.data",
                 lammps_cmd="lammps", db_file=None, parents=None, name="LammpsFW", **kwargs):
        """
        write lammps inputset, run, and store the output.

        Args:
            lammps_input_set (DictLammpsInput): lammps input set
            input_filename (str): input file name
            data_filename (str): data file name
            lammps_cmd (str): command to run lammps.
            db_file (str): path to file specifying db credentials to place output parsing.
            parents ([Fireworks)]: parents of this particular Firework.
            name (str): descriptive name for lammps simulation
            \*\*kwargs: other kwargs that are passed to Firework.__init__.
        """

        t = [
            WriteLammpsFromIOSet(lammps_input_set=lammps_input_set, input_filename=input_filename,
                                 data_filename=data_filename),

            RunLammpsDirect(lammps_cmd=lammps_cmd),

            ToDbTask(drone=LammpsForceFieldDrone(), mmdb="atomate.lammps.database.LammpsCalcDb",
                     db_file=db_file, additional_fields={"task_label": name})
        ]

        super(LammpsFW, self).__init__(t, parents=parents, name=name, **kwargs)


# TODO: implement this
class PackmolFW(Firework):
    pass
