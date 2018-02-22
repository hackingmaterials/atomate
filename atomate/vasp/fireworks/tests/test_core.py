# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.surface import generate_all_slabs
from pymatgen.transformations.standard_transformations import RotationTransformation

from atomate.vasp.fireworks.core import *

__author__ = 'Shyam Dwaraknath'
__email__ = 'shyamd@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
reference_dir = os.path.join(module_dir, "..", "..", "test_files")

DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...


class TestCoreFireworks(unittest.TestCase):
    def setUp(self):
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = Lattice([[3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        self.structure = Structure(lattice, ["Si", "Si"], coords)

    def testOptimizeFW(self):
        self.assertEqual(OptimizeFW(structure=self.structure).name, "Si-structure optimization")

    def testStaticFW(self):
        self.assertEqual(StaticFW(structure=self.structure).name, "Si-static")

        stat_fw = StaticFW(prev_calc_dir="/")
        self.assertEqual(stat_fw.name, "unknown-static")
        self.assertEqual(stat_fw.tasks[0]["calc_dir"], "/")

        opt_fw = OptimizeFW(structure=self.structure)
        stat_fw = StaticFW(parents=opt_fw)
        self.assertEqual(stat_fw.name, "unknown-static")
        self.assertTrue(stat_fw.tasks[0]["calc_loc"])

        self.assertRaises(ValueError, StaticFW)

    def testHSEBSFW(self):
        opt_fw = OptimizeFW(self.structure)
        self.assertEqual(HSEBSFW(structure=self.structure, parents=opt_fw).name, "Si-hse gap")
        self.assertEqual(HSEBSFW(structure=self.structure, parents=opt_fw, mode="line").name, "Si-hse line")

        hse_fw = HSEBSFW(prev_calc_dir="/")
        self.assertEqual(hse_fw.name, "unknown-hse gap")
        self.assertEqual(hse_fw.tasks[0]["calc_dir"], "/")

        hse_fw = HSEBSFW(parents=opt_fw)
        self.assertEqual(hse_fw.name, "unknown-hse gap")
        self.assertTrue(hse_fw.tasks[0]["calc_loc"])

        self.assertRaises(ValueError, HSEBSFW)

    def testNonSCFFW(self):
        opt_fw = OptimizeFW(self.structure)
        self.assertEqual(NonSCFFW(structure=self.structure, parents=opt_fw).name, "Si-nscf uniform")
        self.assertEqual(NonSCFFW(structure=self.structure, parents=opt_fw, mode="line").name, "Si-nscf line")

        nscf_fw = NonSCFFW(prev_calc_dir="/")
        self.assertEqual(nscf_fw.name, "unknown-nscf uniform")
        self.assertEqual(nscf_fw.tasks[0]["calc_dir"], "/")

        nscf_fw = NonSCFFW(parents=opt_fw)
        self.assertEqual(nscf_fw.name, "unknown-nscf uniform")
        self.assertTrue(nscf_fw.tasks[0]["calc_loc"])

        self.assertRaises(ValueError, NonSCFFW)

    def testDFPTFW(self):
        self.assertEqual(DFPTFW(structure=self.structure).name, "Si-static dielectric")

        dfpt_fw = DFPTFW(prev_calc_dir="/")
        self.assertEqual(dfpt_fw.name, "unknown-static dielectric")
        self.assertEqual(dfpt_fw.tasks[0]["calc_dir"], "/")

        opt_fw = OptimizeFW(structure=self.structure)
        dfpt_fw = DFPTFW(parents=opt_fw)
        self.assertEqual(dfpt_fw.name, "unknown-static dielectric")
        self.assertTrue(dfpt_fw.tasks[0]["calc_loc"])

        self.assertRaises(ValueError, DFPTFW)

    def testRamanFW(self):
        opt_fw = OptimizeFW(structure=self.structure)
        self.assertEqual(
            RamanFW(mode=1, displacement=0.01, structure=self.structure, parents=opt_fw).name, "Si-raman_1_0.01")

        raman_fw = RamanFW(mode=1, displacement=0.01, prev_calc_dir="/")
        self.assertEqual(raman_fw.name, "unknown-raman_1_0.01")
        self.assertEqual(raman_fw.tasks[0]["calc_dir"], "/")

        raman_fw = RamanFW(mode=1, displacement=0.01, parents=opt_fw)
        self.assertEqual(raman_fw.name, "unknown-raman_1_0.01")
        self.assertTrue(raman_fw.tasks[0]["calc_loc"])

        self.assertRaises(ValueError, RamanFW, mode=1, displacement=0.01)

    def testSOCFW(self):
        opt_fw = OptimizeFW(self.structure)
        self.assertEqual(SOCFW(magmom=1, structure=self.structure, parents=opt_fw).name, "Si-spin-orbit coupling")

        soc_fw = SOCFW(magmom=1, prev_calc_dir="/")
        self.assertEqual(soc_fw.name, "unknown-spin-orbit coupling")
        self.assertEqual(soc_fw.tasks[0]["calc_dir"], "/")

        soc_fw = SOCFW(magmom=1, parents=opt_fw)
        self.assertEqual(soc_fw.name, "unknown-spin-orbit coupling")
        self.assertTrue(soc_fw.tasks[0]["calc_loc"])

        self.assertRaises(ValueError, SOCFW, magmom=1)

    def testTransmuterFW(self):

        transformations = [RotationTransformation]
        opts = [{"axis": [1, 0, 0], "angle": 0}]

        opt_fw = OptimizeFW(self.structure)
        self.assertEqual(
            TransmuterFW(
                structure=self.structure, transformations=transformations, transformation_params=opts,
                parents=opt_fw).name, "Si-structure transmuter")

    def testBoltztrapFW(self):
        opt_fw = OptimizeFW(self.structure)
        self.assertEqual(BoltztrapFW(structure=self.structure, parents=opt_fw).name, "Si-boltztrap")

        boltz_fw = BoltztrapFW(prev_calc_dir="/")
        self.assertEqual(boltz_fw.name, "unknown-boltztrap")
        self.assertEqual(boltz_fw.tasks[0]["calc_dir"], "/")

        boltz_fw = BoltztrapFW(parents=opt_fw)
        self.assertEqual(boltz_fw.name, "unknown-boltztrap")
        self.assertTrue(boltz_fw.tasks[0]["calc_loc"])

        self.assertRaises(ValueError, BoltztrapFW)

    def testSurfCalcOptimizer(self):

        Fe = Structure.from_spacegroup("Im-3m", Lattice.cubic(2.819), ["Fe"], [(0,0,0)])
        kwargs = {"scratch_dir": ".", "k_product": 50, "db_file": ".",
                  "vasp_cmd": "vasp", "cwd": ".", "mpid": "mp-13"}

        # FW for conventional unit cell calculation
        surface_fw = SurfCalcOptimizer(Fe, structure_type="conventional_unit_cell", **kwargs)
        self.assertEqual(surface_fw.el, "Fe")
        additionals = surface_fw.get_tasks[3]["additional_fields"]
        self.assertEqual(len(additionals.keys()), 5)
        self.assertEqual(additionals["calculation_name"],
                         "Fe_mp-13_conventional_unit_cell_k50")
        self.assertEqual(surface_fw.name, additionals["calculation_name"])
        self.assertEqual(len(surface_fw.get_tasks), 5)

        slabs = generate_all_slabs(Fe, 1, 10, 10, max_normal_search=1,
                                       include_reconstructions=True)
        # FW for oriented unit cell calculation
        slab = slabs[0]
        surface_fw = SurfCalcOptimizer(slab.oriented_unit_cell, structure_type="oriented_unit_cell",
                                       miller_index=slab.miller_index,
                                       scale_factor=slab.scale_factor, **kwargs)

        additionals = surface_fw.get_tasks[3]["additional_fields"]
        self.assertEqual(len(additionals.keys()), 8)
        self.assertEqual(additionals["conventional_spacegroup"],
                         {"symbol": "Im-3m", "number": 229})
        self.assertEqual(additionals["miller_index"], tuple(slab.miller_index))
        self.assertTrue(all(all(t) for t in additionals["scale_factor"] == slab.scale_factor))
        self.assertEqual(additionals["calculation_name"], "Fe_mp-13_bulk_k50_111")
        self.assertEqual(len(surface_fw.get_tasks), 5)

        # FW for reconstructed and unreconstructed slab cell
        slabs = [slabs[0], slabs[-1]]
        for slab in slabs:
            surface_fw = SurfCalcOptimizer(slab, structure_type="slab_cell", vsize=10, ssize=10,
                                           miller_index=slab.miller_index, shift=slab.shift,
                                           scale_factor=slab.scale_factor, ouc=slab.oriented_unit_cell,
                                           reconstruction=slab.reconstruction, **kwargs)

            additionals = surface_fw.get_tasks[4]["additional_fields"]
            self.assertEqual(len(additionals.keys()), 12)
            self.assertEqual(additionals["conventional_spacegroup"],
                             {"symbol": "Im-3m", "number": 229})
            p = list(slab.miller_index)
            p.append(additionals["shift"])
            n = "Fe_mp-13_slab_k50_s10v10_"
            if slab.reconstruction:
                self.assertEqual(additionals["calculation_name"], n+slab.reconstruction)
            else:
                self.assertEqual(additionals["calculation_name"], n+"%s%s%s_shift%s" %tuple(p))
            self.assertEqual(len(surface_fw.get_tasks), 5)



if __name__ == "__main__":
    unittest.main()
