# coding: utf-8

from __future__ import absolute_import, division, print_function, \
    unicode_literals

try:
    from math import gcd
except ImportError:
    from fractions import gcd

import multiprocessing
import os
import json
from fireworks import Firework
from fireworks import ScriptTask, FileWriteTask
from monty.os.path import which
from atomate.utils.utils import get_logger
from atomate.atat.firetasks.parse_outputs import McsqsToDbTask
from atomate.atat.utilities import preprocess_disordered_structure
from fractions import Fraction
from string import Template

__author__ = 'Matthew Horton'
__credits__ = 'Josua Vieten'
__email__ = 'mkhorton@lbl.gov'

logger = get_logger(__name__)
module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

class McsqsFW(Firework):

    def __init__(self, disordered_struct,
                 name="mcsqs",
                 size=None,
                 clusters=None,
                 user_input_settings=None,
                 max_denominator=8,
                 walltime=60,
                 ncores=None,
                 db_file=None,
                 **kwargs):
        """
        Find a best SQS approximation for a given disordered
        structure using the mcsqs code from the ATAT library.

        One unusual feature of this Firework: walltime must
        be specified. This is because there are often cases
        where mcsqs will never terminate (it will keep trying
        structures at random); however, whatever structure
        it finds in this time (the 'bestsqs') is often
        still useful.

        :param disordered_struct: disordered structure
        :param name: name of the Firework
        :param size: ratio of atoms in the output supercell
        to number of atoms in input cell, if None will
        choose the smallest possible supercell
        :param clusters: dict of cluster sizes, e.g.
        {2: 1.5, 3: 1.0} to search for pairs of atoms
        within a radius of 1.5 Å and triplets of atoms
        within a radius of 1.0 Å, if None will try to
        pick sensible defaults, the value of mcsqs is
        that this should converge rapidly with number of
        clusters (aim for around a dozen)
        :param user_input_settings: dict of additional keyword
        value pairs to pass to mcsqs, such as Monte Carlo
        temperature, no validation is performed on this dict
        :param: max_denominator (int):
        Will round occupancies to the nearest rational number,
        with the maximum denominator as specified. Use None not
        to do perform rounding, though not mcsqs might not like this!
        :param: walltime (int): Time to run mcsqs for in minutes,
        1 minute will be deducted from this to run the database
        tasks.
        :param: ncores (int): number of instances of mcsqs
        to run (mcsqs is not parallelized, these are run
        independently of each other), by default will try
        to detect automatically
        :param: db_file:
        :param kwargs: Other kwargs that are passed to Firework.__init__.
        """

        if walltime < 3:
            raise ValueError("Please increase walltime.")

        if disordered_struct.is_ordered:
            raise ValueError("You must input a disordered structure.")

        if db_file is None:
            logger.warn("Database not specified for McsqsFW.")

        # preprocess disordered structure (scales occupancies, normalizes volume)
        disordered_struct = preprocess_disordered_structure(disordered_struct,
                                                            max_denominator=max_denominator)

        if size is None:
            size = self._determine_min_cell(disordered_struct)*2
        else:
            size = self._determine_min_cell(disordered_struct)*size

        if clusters is None:
            # TODO: make cluster determination a bit smarter!
            lattice_param = min(disordered_struct.lattice.abc)
            clusters = {
                -2: lattice_param*2.001,
                -3: lattice_param*2.001,
                -4: lattice_param*1.501
            }
        else:
            for cluster in clusters.keys():
                if cluster not in [2, 3, 4, 5, 6]:
                    raise ValueError("Mcsqs only supports clusters of 2-6 atoms.")

        rndstr = disordered_struct.to(fmt='mcsqs')

        cluster_str = " ".join(["{}={}".format(k, v) for k, v in clusters.items()])
        generate_cluster_cmd = "mcsqs {}".format(cluster_str)

        if user_input_settings:
            user_input_settings_str = " ".join(["{}={}".format(k, v) for k, v
                                            in user_input_settings.items()])
        else:
            user_input_settings_str = ""

        # using same approach as Custodian to detect ncores
        ncores = os.environ.get('NSLOTS') or multiprocessing.cpu_count()
        ncores = int(ncores)

        # mcsqs is not parallelised, and is not expected to finish within the walltime
        # it's a Monte Carlo process that can run indefinitely, so we use timeout
        # to ensure the command finishes
        timeout_cmd = which('timeout') or which('gtimeout')
        if not timeout_cmd:
            raise RuntimeError("Requires timeout to limit walltime.")

        with open(os.path.join(module_dir, 'run_mcsqs_template.sh')) as f:
            template = f.read()

        bash_script = Template(template)
        bash_script = bash_script.substitute({
            'timeout_cmd': timeout_cmd,
            'time': 60*(walltime - 1),  # seconds, minus 60 seconds for DB task
            'ncores': ncores,
            'size': size,
            'settings': user_input_settings_str
        })
        bash_script = str(bash_script)

        # write our input args, so that they can be picked up by the drone
        mcsqs_input_args = json.dumps({
            'clusters': clusters,
            'user_input_settings': user_input_settings,
            'walltime': ncores*(walltime - 1)
        })

        files_to_write = [
            {
                'filename': 'rndstr.in',
                'contents': rndstr
            },
            {
                'filename': 'mcsqs_input_args.json',
                'contents': mcsqs_input_args
            }
        ]

        tasks = [
            FileWriteTask(files_to_write=files_to_write),
            ScriptTask(script=[
                generate_cluster_cmd,
                bash_script
            ], shell_exe='/bin/bash')
        ]

        if db_file:
            tasks.append(McsqsToDbTask(db_file=db_file))

        super(McsqsFW, self).__init__(tasks, name=name, **kwargs)

    @staticmethod
    def _determine_min_cell(disordered_struct):
        """
        Get a minimum cell size for an ordered structure.
        """
        denoms = {}
        for sp_occu in disordered_struct.species_and_occu:
            for sp, occu in sp_occu.items():
                denom = Fraction(occu).limit_denominator(100).denominator
                if sp in denoms:
                    denoms[sp] = max(denom, denoms[sp])
                else:
                    denoms[sp] = denom
        return max(denoms.values())*len(disordered_struct)
