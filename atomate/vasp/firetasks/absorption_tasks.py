from fireworks import FiretaskBase, explicit_serialize
from pymatgen.io.vasp.sets import MPAbsorptionSet


@explicit_serialize
class WriteVaspAbsorptionFromPrev(FiretaskBase):
    """
    Writes input files for an  LOPTICS absorption run. Assumes that output files (WAVECAR) from an
    scf job can be accessed.
    Optional params:
        "prev_calc_dir",
        "mode", either "IPA" or "RPA"
        "reciprocal_density",
        "other_params",
        "potcar_spec"

    """

    optional_params = [
        "prev_calc_dir",
        "structure",
        "mode",
        "copy_wavecar",
        "nbands",
        "nbands_factor",
        "reciprocal_density",
        "nkred",
        "ncores",
        "nedos",
        "potcar_spec",
        "other_params",
    ]

    def run_task(self, fw_spec):
        vis = MPAbsorptionSet.from_prev_calc(
            prev_calc_dir=self.get("prev_calc_dir", "."),
            mode=self.get("mode", "IPA"),
            copy_wavecar=self.get("copy_wavecar", True),
            nbands=self.get("nbands", None),
            nbands_factor=self.get("nbands_factor", 2),
            reciprocal_density=self.get("reciprocal_density", 200),
            nkred=self.get("nkred", None),
            nedos=self.get("nedos", 2001),
            **self.get("other_params", {})
        )
        potcar_spec = self.get("potcar_spec", False)
        vis.write_input(".", potcar_spec=potcar_spec)
