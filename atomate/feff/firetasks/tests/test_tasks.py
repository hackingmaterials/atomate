import os
import unittest
from glob import glob

from pymatgen.core import Structure
from pymatgen.io.feff.inputs import Paths
from pymatgen.io.feff.sets import MPEXAFSSet

from atomate.feff.firetasks.glue_tasks import CopyFeffOutputs
from atomate.feff.firetasks.write_inputs import WriteEXAFSPaths
from atomate.utils.testing import AtomateTest

__author__ = "Kiran Mathew"
__email__ = "kmathew@lbl.gov"

module_dir = os.path.dirname(os.path.abspath(__file__))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")


class TestTasks(AtomateTest):
    def setUp(self):
        super().setUp()
        self.struct = Structure.from_file(
            os.path.join(module_dir, "..", "..", "test_files", "feo_781777.json")
        )

    def test_copy_feff_outputs_task(self):
        t = CopyFeffOutputs(calc_dir=os.path.join(module_dir, "..", "..", "test_files"))
        t.run_task({})
        answer = ["Co2O2.cif", "feff_eels.inp", "feo_781777.json"]
        files = glob("*")
        self.assertEqual(sorted(answer), sorted(files))

    def test_write_paths_task(self):
        exafs = MPEXAFSSet(0, self.struct, edge="K", radius=10)
        t = WriteEXAFSPaths(feff_input_set=exafs, paths=[[249, 0], [85, 0]])
        paths = Paths(exafs.atoms, [[249, 0], [85, 0]])
        paths.write_file("paths_answer.dat")
        t.run_task({})
        with open("paths_answer.dat") as answer, open("paths.dat") as tmp:
            self.assertEqual(answer.readlines(), tmp.readlines())


if __name__ == "__main__":
    unittest.main()
