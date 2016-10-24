# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

from fireworks import Workflow, FileWriteTask
from fireworks.core.firework import Tracker
from fireworks.utilities.fw_utilities import get_slug

from matmethods.utils.utils import get_meta_from_structure, get_fws_and_tasks, update_wf
from matmethods.vasp.firetasks.glue_tasks import CheckStability, CheckBandgap
from matmethods.vasp.firetasks.run_calc import RunVaspCustodian, RunVaspDirect, RunVaspFake
from matmethods.vasp.firetasks.write_inputs import ModifyIncar
from matmethods.vasp.config import ADD_NAMEFILE, SCRATCH_DIR, ADD_MODIFY_INCAR

__author__ = 'Anubhav Jain, Kiran Mathew'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov'

import warnings

warnings.warn("vasp_config renamed to config.")

from .powerups import *