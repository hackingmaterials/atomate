# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
import json

from fireworks import FWorker
from fireworks.core.rocket_launcher import rapidfire
from atomate.utils.testing import AtomateTest
from pymatgen.io.qchem.outputs import QCOutput
from atomate.qchem.workflows.base.ion_placement import get_ion_placement_wf
from atomate.qchem.database import QChemCalcDb
from pymatgen.core import Molecule
try:
    from unittest.mock import patch, MagicMock
except ImportError:
    from mock import patch, MagicMock

__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "11/15/18"


module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
test_file_dir = os.path.join(module_dir, "..", "..","test_files","ion_placer_files")

class TestIonPlacement(AtomateTest):
    def test_Ion_Placement(self):
        ec_out = QCOutput(os.path.join(test_file_dir,"EC-12.qout")).data
        mol = ec_out['initial_molecule']
        mol.add_site_property("charge",ec_out['Mulliken'][0][::,0])
        mol.add_site_property("spin",ec_out['Mulliken'][0][::,1])
        test_positions = [[-2.5493143265,0.5101339356,-1.7108711165],[-0.6615990397,2.3695796927,0.9279735702],[-0.3750215537,-0.0304037817,-1.7194196010],[0.7930586053,0.0296565175,1.9128737508],[-1.9492306263,-3.1583707935,-0.1981875209],[0.5401522606,1.7637013727,-0.6447137855],[1.8492267629,0.1903840521,-0.2156747400],[-1.6463297659,1.2830427999,2.2135363340],[-2.5586859415,2.7187118030,0.5463830493],[-1.7466289410,-2.2770280050,-1.8155787349],[1.1548402487,-1.6224314280,-0.1450955593],[-1.3435798028,2.4808810589,-0.8023440594],[-1.2021394971,-1.4831508185,1.4310011072],[-5.2837371547,-0.0455308548,-0.1321178451]]
        ref_dirs = {}
        for ii in range(14):
            ref_dirs["ion_pos_"+str(ii)] = os.path.join(test_file_dir,"launcher"+str(ii))
        wf = get_ion_placement_wf(
            molecule=mol,
            ion="Li",
            pcm_dielectric=30.0,
            do_optimization=False,
            qchem_input_params={"overwrite_inputs": {"rem": {"XC_GRID": "3"}}},
            test_positions=test_positions,
            ref_dirs=ref_dirs)
        self.lp.add_wf(wf)
        rapidfire(self.lp,fworker=FWorker(env={"max_cores": 32, "db_file": os.path.join(db_dir, "db.json")}), pdb_on_exception=True)
        mmdb = QChemCalcDb.from_db_file(os.path.join(db_dir, "db.json"), admin=True)
        target_entries = list(mmdb.collection.find({"task_label": "gather_geoms"}))
        self.assertEqual(target_entries[0]["sorted_data"][0]["energy"],-349.947920236903)
        mmdb.reset()

if __name__ == "__main__":
    unittest.main()
