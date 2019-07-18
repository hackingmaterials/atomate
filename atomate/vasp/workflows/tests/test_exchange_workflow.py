# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

from monty.os.path import which

from atomate.vasp.workflows.base.exchange import ExchangeWF
from atomate.vasp.firetasks.parse_outputs import (
    MagneticDeformationToDB,
    MagneticOrderingsToDB,
)
from atomate.utils.testing import AtomateTest, DB_DIR

from json import load
from pymatgen import Structure

__author__ = "Nathan C. Frey"
__email__ = "ncfrey@lbl.gov"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
ref_dir = os.path.join(module_dir, "..", "..", "test_files", "exchange_wf")

enum_cmd = which("enum.x") or which("multienum.x")
VAMPEXE = which("vampire-serial")

class TestExchangeWF(AtomateTest):

    @unittest.skipIf(not enumlib_present, "enumlib not present")
    def test_workflow(self):

        # Get WF with static orderings, Heisenber model mapping, and Monte Carlo calculation of critical temp.

        pass
        
