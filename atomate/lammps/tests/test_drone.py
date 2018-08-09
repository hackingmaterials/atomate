# coding: utf-8
from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

# from pymatgen.io.lammps.sets import LammpsInputSet
# from pymatgen.io.lammps.output import LammpsLog

from atomate.utils.testing import AtomateTest
from atomate.lammps.drones import LammpsDrone

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "common", "test_files")

LAMMPS_CMD = None  # "lmp_serial"


@unittest.skip("Atomate Lammps does not work after pymatgen release v2018.5.22")
class TestLammpsDrone(AtomateTest):

    def setUp(self):
        super(TestLammpsDrone, self).setUp()
        self.drone = LammpsDrone()
        self.calc_dir = os.path.join(module_dir, "test_files")
        self.input_file = os.path.join(self.calc_dir, "lammps.in")
        self.data_file = os.path.join(self.calc_dir, "lammps.data")
        self.log_file = os.path.join(self.calc_dir, "lammps.log")

    def test_assimilate(self):
        doc = self.drone.assimilate(self.calc_dir, input_filename="lammps.in",
                                    log_filename="lammps.log",
                                    data_filename="lammps.data")
        lmps_input_set = LammpsInputSet.from_file("lammps", self.input_file,
                                                  {}, lammps_data=self.data_file,
                                                  data_filename="lammps.data")
        # no dump file ==> output is just the log file
        lmps_output = LammpsLog(self.log_file)
        self.assertDictEqual(doc["input"], lmps_input_set.as_dict())
        self.assertDictEqual(doc["output"]["log"], lmps_output.as_dict())

        enthalpy = [1906.1958, 1220.2265, 596.51973, 465.01619, 148.91822, 26.160144,
                    319.27146, 141.35729, 299.04503, 271.19625, 145.4361]
        self.assertEqual(lmps_output.as_dict()['thermo_data']['enthalpy'], enthalpy)


if __name__ == "__main__":
    unittest.main()
