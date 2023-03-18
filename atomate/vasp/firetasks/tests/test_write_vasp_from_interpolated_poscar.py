import os
import unittest

from fireworks import Firework, Workflow
from fireworks.core.rocket_launcher import rapidfire
from pymatgen.core.structure import Structure

from atomate.common.firetasks.glue_tasks import PassCalcLocs, get_calc_loc
from atomate.utils.testing import AtomateTest
from atomate.vasp.firetasks.glue_tasks import CopyVaspOutputs
from atomate.vasp.firetasks.write_inputs import WriteVaspFromIOSetFromInterpolatedPOSCAR

__author__ = "Tess Smidt"
__email__ = "blondegeek@gmail.com"

module_dir = os.path.dirname(os.path.abspath(__file__))

DEBUG_MODE = False


class TestWriteVaspFromInterpolatedPOSCAR(AtomateTest):
    @classmethod
    def setUpClass(cls):
        cls.static_outdir = os.path.join(
            module_dir, "..", "..", "test_files", "Si_static", "outputs"
        )
        cls.opt_outdir = os.path.join(
            module_dir, "..", "..", "test_files", "Si_structure_optimization", "outputs"
        )

    def test_writevaspfrominterpolatedposcar(self):
        nimages = 5
        this_image = 1
        autosort_tol = 0.5

        fw1 = Firework(
            [
                CopyVaspOutputs(
                    calc_dir=self.static_outdir,
                    contcar_to_poscar=False,
                    additional_files=["CONTCAR"],
                ),
                PassCalcLocs(name="fw1"),
            ],
            name="fw1",
        )

        fw2 = Firework(
            [
                CopyVaspOutputs(
                    calc_dir=self.opt_outdir,
                    contcar_to_poscar=False,
                    additional_files=["CONTCAR"],
                ),
                PassCalcLocs(name="fw2"),
            ],
            name="fw2",
        )

        fw3 = Firework(
            [
                WriteVaspFromIOSetFromInterpolatedPOSCAR(
                    start="fw1",
                    end="fw2",
                    this_image=this_image,
                    nimages=nimages,
                    autosort_tol=autosort_tol,
                    vasp_input_set="MPStaticSet",
                ),
                PassCalcLocs(name="fw3"),
            ],
            name="fw3",
            parents=[fw1, fw2],
        )

        fw4 = Firework([PassCalcLocs(name="fw4")], name="fw4", parents=fw3)

        wf = Workflow([fw1, fw2, fw3, fw4])
        self.lp.add_wf(wf)
        rapidfire(self.lp)

        fw4 = self.lp.get_fw_by_id(self.lp.get_fw_ids({"name": "fw4"})[0])

        calc_locs = fw4.spec["calc_locs"]

        print(get_calc_loc("fw3", calc_locs)["path"])

        # Check existence of structure files.
        self.assertTrue(
            os.path.exists(get_calc_loc("fw3", calc_locs)["path"] + "/POSCAR")
        )
        self.assertTrue(
            os.path.exists(
                get_calc_loc("fw3", calc_locs)["path"] + "/interpolate/CONTCAR_0"
            )
        )
        self.assertTrue(
            os.path.exists(
                get_calc_loc("fw3", calc_locs)["path"] + "/interpolate/CONTCAR_1"
            )
        )

        self.assertTrue(
            os.path.exists(get_calc_loc("fw3", calc_locs)["path"] + "/INCAR")
        )
        self.assertTrue(
            os.path.exists(get_calc_loc("fw3", calc_locs)["path"] + "/KPOINTS")
        )
        self.assertTrue(
            os.path.exists(get_calc_loc("fw3", calc_locs)["path"] + "/POTCAR")
        )

        # Check interpolation.
        struct_start = Structure.from_file(
            get_calc_loc("fw3", calc_locs)["path"] + "/interpolate/CONTCAR_0"
        )
        struct_end = Structure.from_file(
            get_calc_loc("fw3", calc_locs)["path"] + "/interpolate/CONTCAR_1"
        )
        struct_inter = Structure.from_file(
            get_calc_loc("fw3", calc_locs)["path"] + "/POSCAR"
        )

        structs = struct_start.interpolate(
            struct_end, nimages, interpolate_lattices=True, autosort_tol=autosort_tol
        )

        # Check x of 1st site.
        self.assertAlmostEqual(
            structs[this_image][1].coords[0], struct_inter[1].coords[0]
        )
        # Check c lattice parameter
        self.assertAlmostEqual(
            structs[this_image].lattice.abc[0], struct_inter.lattice.abc[0]
        )


if __name__ == "__main__":
    unittest.main()
