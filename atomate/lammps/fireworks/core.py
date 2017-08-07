# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
Defines fireworks to be incorporated into workflows.
"""

from fireworks import Firework

from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.lammps.firetasks.run_calc import RunLammpsDirect, RunPackmol
from atomate.lammps.firetasks.parse_outputs import LammpsToDB
from atomate.lammps.firetasks.write_inputs import WriteFromIOSet, WriteFromForceFieldAndTopology
from atomate.lammps.firetasks.glue_tasks import CopyPackmolOutputs

__author__ = "Brandon Wood, Kiran Mathew"
__email__ = "b.wood@berkeley.edu"


class LammpsFW(Firework):

    def __init__(self, lammps_input_set, input_filename="lammps.in", data_filename="lammps.data",
                 lammps_cmd="lammps", db_file=None, parents=None, name="LammpsFW",
                 log_filename="log.lammps", dump_filename=None, **kwargs):
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

        tasks = [
            WriteFromIOSet(lammps_input_set=lammps_input_set, input_filename=input_filename,
                           data_filename=data_filename),

            RunLammpsDirect(lammps_cmd=lammps_cmd, input_filename=input_filename),

            LammpsToDB(input_filename=input_filename, data_filename=data_filename,
                       log_filename=log_filename, dump_filename=dump_filename,
                       db_file=db_file, additional_fields={"task_label": name})
        ]

        super(LammpsFW, self).__init__(tasks, parents=parents, name=name, **kwargs)


class LammpsForceFieldFW(Firework):

    def __init__(self, input_file, final_molecule_path, molecules, mols_number, forcefield,
                 user_settings=None, site_property=None, input_filename="lammps.in",
                 data_filename="lammps.data", lammps_cmd="lammps", db_file=None, parents=None,
                 log_filename="log.lammps", dump_filename=None, name="LammpsFFFW", **kwargs):

        user_settings = user_settings or {}

        tasks = [
            CopyPackmolOutputs(calc_loc=True),

            WriteFromForceFieldAndTopology(input_file=input_file, final_molecule_path=final_molecule_path,
                                           molecules=molecules, mols_number=mols_number,
                                           forcefield=forcefield, input_filename=input_filename,
                                           user_settings=user_settings, site_property=site_property),

            RunLammpsDirect(lammps_cmd=lammps_cmd, input_filename=input_filename),

            LammpsToDB(input_filename=input_filename, data_filename=data_filename,
                       log_filename=log_filename, dump_filename=dump_filename,
                       db_file=db_file, additional_fields={"task_label": name})
        ]

        super(LammpsForceFieldFW, self).__init__(tasks, parents=parents, name=name, **kwargs)


class PackmolFW(Firework):

    def __init__(self, molecules, packing_config, tolerance=2.0, filetype="xyz", control_params=None,
                 output_file="packed.xyz",  copy_to_current_on_exit=False, site_property=None,
                 parents=None, name="PackmolFW", **kwargs):
        """

        Args:
            molecules:
            packing_config:
            tolerance:
            filetype:
            control_params:
            output_file:
            parents:
            name:
            **kwargs:
        """
        control_params = control_params or {'maxit': 20, 'nloop': 600}

        tasks = [
            RunPackmol(molecules=molecules, packing_config=packing_config, tolerance=tolerance,
                       filetype=filetype, control_params=control_params,  output_file=output_file,
                       copy_to_current_on_exit=copy_to_current_on_exit, site_property=site_property),

            PassCalcLocs(name=name)

             ]

        super(PackmolFW, self).__init__(tasks, parents=parents, name=name, **kwargs)
