# coding: utf-8

from __future__ import absolute_import, division, print_function, \
    unicode_literals

try:
    from math import gcd
except ImportError:
    from fractions import gcd

import warnings
import multiprocessing
import os
import numpy as np
from fireworks import Firework
from fireworks import ScriptTask
from pymatgen import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher, OccupancyComparator
from monty.serialization import dumpfn
from pymatgen.transformations.advanced_transformations import DiscretizeOccupanciesTransformation
from fractions import Fraction
from pymatgen.transformations.advanced_transformations import EnumerateStructureTransformation
from pymatgen.transformations.standard_transformations import OrderDisorderedStructureTransformation


class McsqsFW(Firework):

    def __init__(self, disordered_struct,
                 path,
                 name="mcsqs",
                 size=None,
                 sizemultipl=4,
                 max_den=8,
                 tolerance=0.5,
                 clusters=None,
                 user_input_settings=None,
                 occu_tol=8,
                 walltime=2.5,
                 ncores=None,
                 **kwargs):
        """
        Find a best SQS approximation for a given disordered
        structure using the mcsqs code from the ATAT library.

        :param disordered_struct: disordered structure
        :param name: name of the Firework
        :param size: number of atoms in the output supercell,
        if None will choose the smallest possible supercell
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
        :param: occu_tol (int): As in pymatgen's
        EnumerateStructureTransformation.
        Will round occupancies to the nearest rational number,
        with the maximum denominator equal to occu_tol.
        :param: walltime (int): Time to run mcsqs for in minutes,
        2 minutes will be deducted from this to run the database
        tasks.
        :param: ncores (int): number of instances of mcsqs
        to run (mcsqs is not parallelized, these are run
        independently of each other), by default will try
        to detect automatically
        :param kwargs: Other kwargs that are passed to Firework.__init__.
        """

        if disordered_struct.is_ordered:
            raise ValueError("You must input a disordered structure.")

        if occu_tol:
            disordered_struct = self._discretize_cell(disordered_struct,max_denominator=occu_tol)

        if size is None:
            size = self._determine_min_cell(disordered_struct)   #this would be a more elegant way to determine the size, but it seems to take forever

        if clusters is None:
            # TODO: make cluster determination a bit smarter!
            lattice_param = min(disordered_struct.lattice.abc)
            clusters = {
                -2: lattice_param*1.5,
                -3: lattice_param*1.5,
                -4: lattice_param
            }
        else:
            for cluster in clusters.keys():
                if cluster not in [2, 3, 4, 5, 6]:
                    raise ValueError("Mcsqs only supports clusters of 2-6 atoms.")


        disordered_struct.to(filename='rndstr2.in')   

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

        run_mcsqs_cmd = """
for (( id=0 ; id<{} ; id ++ ))

do
    timeout {}m mcsqs -n {} {}  -ip=$id & 
done
""".format(ncores, walltime-2, size, user_input_settings_str)

        # find the best SQS among all that have been run in parallel (using mcsqs' own method; lowest objective function)
        get_bestsqs_cmd = "sleep " + str(walltime) + "m; mcsqs -best" 

        # copy files from the original directory into the working directory
        copyfiles = "cp -r -f " + str(path) + "* ."

        # copy resulting files from the working directory back into the originial directory
        copyfiles1 = "cp -f *.out "  + str(path) + " && cp -f *.log "  + str(path) + " && cp -f FW.json "  + str(path) + " && cp -f clusters.out "  + str(path) + " && cp -f mcsqs_input_args.json "  + str(path) + " && cp -f mcsqs_version.txt "  + str(path)

        # write the mcsqs version to a file for provenance
        write_version_cmd = "mcsqs -v &> mcsqs_version1.txt; head -n 1 mcsqs_version1.txt > mcsqs_version.txt; rm -f mcsqs_version1.txt"

        # write our input args, so that they can be picked up by the drone
        dumpfn({
            'clusters': clusters,
            'user_input_settings': user_input_settings,
            'walltime': ncores*(walltime-2)
        }, "mcsqs_input_args.json")

        # do not generate clusters if they have already been generated
        try: 
            filen = open(str(path) + "clusters.out")
            cl = 1
        except:
            cl = 0

        if cl == 0:
            t = [
                ScriptTask(script=[
                    copyfiles,
                    generate_cluster_cmd,
                    run_mcsqs_cmd,
                    get_bestsqs_cmd,
                    write_version_cmd,
                    copyfiles1
                ], shell_exe='/bin/bash')
            ]

        else:

            t = [
                ScriptTask(script=[
                    copyfiles,
                    run_mcsqs_cmd,
                    get_bestsqs_cmd,
                    write_version_cmd,
                    copyfiles1
                ], shell_exe='/bin/bash')
            ]

        # TODO: add tracker(s)

        super(McsqsFW, self).__init__(t, name=name, **kwargs)

    @staticmethod
    def _discretize_cell(disordered_struct,max_denominator=8):
        """
        Remove partial occupancies.
        """
        trans = DiscretizeOccupanciesTransformation(max_denominator)
        disc_struct = trans.apply_transformation(disordered_struct)
        return disc_struct

    @staticmethod
    def _determine_min_cell(disordered_struct):
        """
        Enumerate superstructures and determine the smallest superstructure with discrete occupancies.
        Seems to be very slow.
        """
        return len(disordered_struct)*8
               

# TODO: add a duplicate checker
# when we want to create a firework, check the database first
# if bestsqs present, use that, if not run McsqsFW()

# to do this, will query database for all structures with same space group
# then for every disordered structure in returned query, run structurematcher
# if we get a match, use the bestsqs from that doc
# sm = StructureMatcher(comparator=OccupancyComparator())
