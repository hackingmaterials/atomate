# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines tasks that support running FEFF.
"""

import subprocess

from fireworks import explicit_serialize, FireTaskBase
import os
import sys
import subprocess

from matmethods.utils.utils import env_chk, get_logger

__author__ = 'Kiran Mathew'
__email__ = "kmathew@lbl.gov"
logger = get_logger(__name__)



SUBMIT_FNAME = "submit_script"


TEMPLATE = """#!/bin/bash
#PBS -q {queue}
#PBS -N {name}
#PBS -l nodes=1:ppn={nproc}:{proc}
#PBS -l walltime={walltime}
#PBS -o $PBS_JOBID.log
#PBS -e {name}.err
#PBS -V
#PBS -M {user}@ucsd.edu
#PBS -m {verbosity}
#PBS -A ong-group
#PBS -d {dir}


mpirun -np {nproc} --hostfile $PBS_NODEFILE {FEFF_HOME}/rdinp
mpirun -np {nproc} --hostfile $PBS_NODEFILE {FEFF_HOME}/atomic
mpirun -np {nproc} --hostfile $PBS_NODEFILE {FEFF_HOME}/dmdw
mpirun -np {nproc} --hostfile $PBS_NODEFILE {FEFF_HOME}/pot
mpirun -np {nproc} --hostfile $PBS_NODEFILE {FEFF_HOME}/ldos
mpirun -np {nproc} --hostfile $PBS_NODEFILE {FEFF_HOME}/screen
mpirun -np {nproc} --hostfile $PBS_NODEFILE {FEFF_HOME}/opconsat
mpirun -np {nproc} --hostfile $PBS_NODEFILE {FEFF_HOME}/xsph
mpirun -np {nproc} --hostfile $PBS_NODEFILE {FEFF_HOME}/fms
mpirun -np {nproc} --hostfile $PBS_NODEFILE {FEFF_HOME}/mkgtr
mpirun -np {nproc} --hostfile $PBS_NODEFILE {FEFF_HOME}/path
mpirun -np {nproc} --hostfile $PBS_NODEFILE {FEFF_HOME}/genfmt
mpirun -np {nproc} --hostfile $PBS_NODEFILE {FEFF_HOME}/ff2x
mpirun -np {nproc} --hostfile $PBS_NODEFILE {FEFF_HOME}/sfconv
mpirun -np {nproc} --hostfile $PBS_NODEFILE {FEFF_HOME}/compton
mpirun -np {nproc} --hostfile $PBS_NODEFILE {FEFF_HOME}/eels


"""

walltime_settings = {
        "home": (24, 240),
        "hotel": (24, 168),
        "condo": (8, 8),
        "glean": (24, 240),
    }

queue_submission_cmd = "qsub"
pjoin = os.path.join


def proc_submit_script(d, FEFF_HOME, name,queue,walltime,proc,verbosity):

    p= {
        "queue":queue,
        "name":name,
        "user":os.environ["USER"],
        "proc":proc,
        "verbosity":verbosity,
        "dir":os.path.abspath(d),
        "FEFF_HOME":FEFF_HOME
    }

    if proc == "haswell":
        p["nproc"] = 24
    elif proc =="sandy":
        p["nproc"] = 16

    p["walltime"]="%d:00:00"%walltime

    with open(os.path.join(d,SUBMIT_FNAME),"w") as f:
        f.write(TEMPLATE.format(**p))




@explicit_serialize
class RunFeffTscc(FireTaskBase):
    """
    Run FEFF directly on TSCC (no custodian)

    Need to have environment variable FEFF_HOME set properly
    """
    optional_params = ['walltime', 'name', 'verbosity', 'queue', 'proc']

    def run_task(self, fw_spec):

        FEFF_HOME = os.environ.get("FEFF_HOME")
        walltime = 4
        name = "feff_job"
        queue = "glean"
        proc = "haswell"
        verbosity = "a"
        d = os.getcwd()


        if "walltime" in self:
            walltime = self["walltime"]

        if "name" in self:
            name = self["name"]

        if "verbosity" in self:
            verbosity = self["verbosity"]

        if "queue" in self:
            queue = self["queue"]

        if "proc" in self:
            proc = self["proc"]

        proc_submit_script(d,FEFF_HOME,name,queue,walltime,proc,verbosity)
        return_code = subprocess.call([queue_submission_cmd, SUBMIT_FNAME])
        logger.info("FEFF finished running with returncode: {}".format(return_code))







