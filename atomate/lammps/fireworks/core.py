# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
Defines fireworks to be incorporated into workflows.
"""

from fireworks import Firework

from atomate.lammps.firetasks.run_calc import RunLammpsDirect, RunPackmol
from atomate.lammps.firetasks.parse_outputs import LammpsToDB
from atomate.lammps.firetasks.write_inputs import WriteLammpsFromIOSet

__author__ = "Brandon Wood, Kiran Mathew"
__email__ = "b.wood@berkeley.edu"


class LammpsFW(Firework):

    def __init__(self, lammps_input_set, input_filename="lammps.in", data_filename="lammps.data",
                 lammps_cmd="lammps", db_file=None, parents=None, name="LammpsFW",
                 log_filename="lammps.log", dump_filename=None, **kwargs):
        """
        write lammps inputset, run, and store the output.

        Args:
            lammps_input_set (DictLammpsInput): lammps input set
            input_filename (str): input file name
            data_filename (str): data file name
            lammps_cmd (str): command to run lammps (skip the input file).
                e.g. 'mpirun -n 8 lmp_mpi'
            db_file (str): path to file specifying db credentials to place output parsing.
            parents ([Fireworks)]: parents of this particular Firework.
            name (str): descriptive name for lammps simulation
            \*\*kwargs: other kwargs that are passed to Firework.__init__.
        """

        t = [
            WriteLammpsFromIOSet(lammps_input_set=lammps_input_set, input_filename=input_filename,
                                 data_filename=data_filename),

            RunLammpsDirect(lammps_cmd=lammps_cmd, input_filename=input_filename),

            LammpsToDB(input_filename=input_filename, data_filename=data_filename,
                       log_filename=log_filename, dump_filename=dump_filename,
                       db_file=db_file, additional_fields={"task_label": name})
        ]

        super(LammpsFW, self).__init__(t, parents=parents, name=name, **kwargs)


class PackmolFW(Firework):

    def __init__(self, molecules, packing_config, tolerance=2.0, filetype="xyz", control_params=None,
                 output_file="packed.xyz", parents=None, name="PackmolFW", **kwargs):
        control_params = control_params or {'maxit': 20, 'nloop': 600}
        t = [
            RunPackmol(molecules=molecules, packing_config=packing_config, tolerance=tolerance,
                       filetype=filetype, control_params=control_params,  output_file=output_file)
             ]
        super(PackmolFW, self).__init__(t, parents=parents, name=name, **kwargs)
