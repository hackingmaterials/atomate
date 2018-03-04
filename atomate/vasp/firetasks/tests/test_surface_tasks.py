# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

from atomate.vasp.firetasks.surface_tasks import FacetFWsGeneratorTask

from pymatgen.io.vasp import Poscar
from pymatgen import Structure, Lattice
from pymatgen.core.surface import SlabGenerator

__author__ = 'Richard Tran'
__email__ = 'rit001@eng.ucsd.edu'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class TestFacetFWsGeneratorTask(unittest.TestCase):
    def setUp(self):

        self.Fe = Structure(Lattice.cubic(2.819), ["Fe", "Fe"],
                            [(0,0,0), (0.5,0.5,0.5)])
        self.Co = Structure(Lattice.hexagonal(2.5, 4.07), ["Co", "Co"],
                            [[1/3, 2/3, 1/4], [2/3, 1/3, 3/4]])
        os.chdir(os.path.join(module_dir, "../../test_files/surface_wf"))

    def test_get_oriented_ucell_fw(self):

        # Test for conventional unit cell
        p = Poscar(self.Fe)
        p.write_file("CONTCAR.relax2.gz")

        facettask = FacetFWsGeneratorTask(structure_type="conventional_unit_cell",
                                          scratch_dir=".", k_product=50, db_file=".",
                                          vasp_cmd="vasp", max_index=1, naming_tag="mp-13")

        fwaction = facettask.run_task({})

        self.assertEqual(len(fwaction.additions), 4)
        for fw in fwaction.additions:
            self.assertEqual(len(fw.tasks), 5)
        all_names = tuple(fw.name for fw in fwaction.additions)
        names = ["Fe_mp-13_bulk_rec_k50_bcc_100_zigzag_rt2xrt2",
                 "Fe_mp-13_bulk_k50_111", "Fe_mp-13_bulk_k50_100",
                 "Fe_mp-13_bulk_k50_110"]
        for n in all_names:
            self.assertTrue(n in names)


    def test_get_slab_fw(self):

        # Test for oriented unit cell
        oriented_ucell = SlabGenerator(self.Co, (1, 0, 2), 10, 10,
                            max_normal_search=2).get_slab().oriented_unit_cell
        p = Poscar(oriented_ucell)
        p.write_file("CONTCAR.relax2.gz")

        facettask = FacetFWsGeneratorTask(structure_type="oriented_unit_cell",
                                          scratch_dir=".", k_product=50,
                                          db_file=".", vasp_cmd="vasp",
                                          miller_index=(1, 0, 2), naming_tag="mp-54")

        fwaction = facettask.run_task({})

        self.assertEqual(len(fwaction.additions), 2)
        for fw in fwaction.additions:
            self.assertEqual(len(fw.tasks), 5)
        all_names = tuple(fw.name for fw in fwaction.additions)
        for n in all_names:
            self.assertTrue("Co_mp-54_slab_k50_s10v10_102" in n)

    def tearDown(self):
        os.remove("CONTCAR.relax2.gz")


if __name__ == '__main__':
    unittest.main()
