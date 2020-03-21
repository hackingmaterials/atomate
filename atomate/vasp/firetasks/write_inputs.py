"""
This module defines tasks for writing vasp input sets.
"""

import os
from importlib import import_module

import numpy as np

from monty.serialization import dumpfn

from fireworks import FiretaskBase, explicit_serialize
from fireworks.utilities.dict_mods import apply_mod

from pymatgen.core.structure import Structure
from pymatgen.alchemy.materials import TransformedStructure
from pymatgen.alchemy.transmuters import StandardTransmuter
from pymatgen.io.vasp import (
    Incar,
    Poscar,
    Potcar,
    PotcarSingle,
    Vasprun,
    Kpoints
)
from pymatgen.io.vasp.sets import (
    MPStaticSet,
    MPNonSCFSet,
    MPSOCSet,
    MPHSEBSSet,
    MPNMRSet,
    MPScanRelaxSet
)

from atomate.utils.utils import env_chk, load_class
from atomate.vasp.firetasks.glue_tasks import GetInterpolatedPOSCAR

__author__ = "Anubhav Jain, Shyue Ping Ong, Kiran Mathew, Alex Ganose"
__email__ = "ajain@lbl.gov"


@explicit_serialize
class WriteVaspFromIOSet(FiretaskBase):
    """
    Create VASP input files using implementations of pymatgen's
    AbstractVaspInputSet. An input set can be provided as an object or as a
    String/parameter combo.

    Required params:
        structure (Structure): structure
        vasp_input_set (AbstractVaspInputSet or str): Either a VaspInputSet
            object or a string name for the VASP input set (e.g., "MPRelaxSet").

    Optional params:
        vasp_input_params (dict): When using a string name for VASP input set,
            use this as a dict to specify kwargs for instantiating the input set
            parameters. For example, if you want to change the
            user_incar_settings, you should provide:
            {"user_incar_settings": ...}. This setting is ignored if you provide
            the full object representation of a VaspInputSet rather than a
            String.
        potcar_spec (bool): Instead of writing the POTCAR, write a
            "POTCAR.spec". This is intended to allow testing of workflows
            without requiring pseudo-potentials to be installed on the system.
    """

    required_params = ["structure", "vasp_input_set"]
    optional_params = ["vasp_input_params", "potcar_spec"]

    def run_task(self, fw_spec):
        # if a full VaspInputSet object was provided
        if hasattr(self["vasp_input_set"], "write_input"):
            vis = self["vasp_input_set"]

        # if VaspInputSet String + parameters was provided
        else:
            vis_cls = load_class(
                "pymatgen.io.vasp.sets", self["vasp_input_set"]
            )
            vis = vis_cls(
                self["structure"], **self.get("vasp_input_params", {})
            )

        potcar_spec = self.get("potcar_spec", False)
        vis.write_input(".", potcar_spec=potcar_spec)


@explicit_serialize
class WriteVaspFromIOSetFromInterpolatedPOSCAR(GetInterpolatedPOSCAR):
    """
    Grabs CONTCARS from two previous calculations to create interpolated
    structure. Create VASP input files using implementations of pymatgen's
    AbstractVaspInputSet. An input set can be provided as String/parameter
    combo.

    Required params:
        start (str): name of fw for start of interpolation.
        end (str): name of fw for end of interpolation.
        this_image (int): which interpolation this is.
        nimages (int) : number of interpolations.
        autosort_tol (float): a distance tolerance in angstrom in which
          to automatically sort end_structure to match to the closest
          points in this particular structure.
        vasp_input_set (str): a string name for the VASP input set (e.g.,
            "MPRelaxSet").

    Optional params:
        vasp_input_params (dict): When using a string name for VASP input set,
            use this as a dict to specify kwargs for instantiating the input set
            parameters. For example, if you want to change the
            user_incar_settings, you should provide:
            {"user_incar_settings": ...}. This setting is ignored if you provide
            the full object representation of a VaspInputSet rather than a
            String.
        potcar_spec (bool): Instead of writing the POTCAR, write a
            "POTCAR.spec". This is intended to allow testing of workflows
            without requiring pseudo-potentials to be installed on the system.
    """

    # First, we make a fresh copy of the required_params before modifying.
    required_params = [
        "start",
        "end",
        "this_image",
        "nimages",
        "vasp_input_set",
    ]
    optional_params = ["vasp_input_params", "autosort_tol", "potcar_spec"]

    def run_task(self, fw_spec):
        # Get interpolated structure.
        structure = GetInterpolatedPOSCAR.interpolate_poscar(self, fw_spec)

        # Assumes VaspInputSet String + parameters was provided
        vis_cls = load_class("pymatgen.io.vasp.sets", self["vasp_input_set"])
        vis = vis_cls(structure, **self.get("vasp_input_params", {}))
        potcar_spec = self.get("potcar_spec", False)
        vis.write_input(".", potcar_spec=potcar_spec)


@explicit_serialize
class WriteVaspFromPMGObjects(FiretaskBase):
    """
    Write VASP files using pymatgen objects.

    Note, that although this firetask has no required params, it is
    non-functional unless at least one optional param is set.

    Optional params:
        incar (Incar): pymatgen Incar object
        poscar (Poscar): pymatgen Poscar object
        kpoints (Kpoints): pymatgen Kpoints object
        potcar (Potcar): pymatgen Potcar object
    """

    optional_params = ["incar", "poscar", "kpoints", "potcar"]

    def run_task(self, fw_spec):
        if "incar" in self:
            self["incar"].write_file("INCAR")
        if "poscar" in self:
            self["poscar"].write_file("POSCAR")
        if "kpoints" in self:
            self["kpoints"].write_file("KPOINTS")
        if "potcar" in self:
            self["potcar"].write_file("POTCAR")


@explicit_serialize
class ModifyIncar(FiretaskBase):
    """
    Modify an INCAR file.

    Required params:
        (none)

    Optional params:
        incar_update (dict): overwrite Incar dict key. Supports env_chk.
        incar_multiply ([{<str>:<float>}]) - multiply Incar key by a constant
            factor. Supports env_chk.
        incar_dictmod ([{}]): use DictMod language to change Incar.
            Supports env_chk.
        input_filename (str): Input filename (if not "INCAR")
        output_filename (str): Output filename (if not "INCAR")
    """

    optional_params = [
        "incar_update",
        "incar_multiply",
        "incar_dictmod",
        "input_filename",
        "output_filename",
    ]

    def run_task(self, fw_spec):

        incar_name = self.get("input_filename", "INCAR")
        incar = Incar.from_file(incar_name)

        incar_update = env_chk(self.get("incar_update"), fw_spec)
        incar_multiply = env_chk(self.get("incar_multiply"), fw_spec)
        incar_dictmod = env_chk(self.get("incar_dictmod"), fw_spec)

        if incar_update:
            incar.update(incar_update)

        if incar_multiply:
            for k in incar_multiply:
                if hasattr(incar[k], "__iter__"):  # is list-like
                    incar[k] = list(np.multiply(incar[k], incar_multiply[k]))
                else:
                    incar[k] = incar[k] * incar_multiply[k]

        if incar_dictmod:
            apply_mod(incar_dictmod, incar)

        incar.write_file(self.get("output_filename", "INCAR"))


@explicit_serialize
class ModifyKpoints(FiretaskBase):
    """
    Modify an KPOINTS file.

    Required params:
        (none)

    Optional params:
        kpoints_update (dict): overwrite Kpoint dict key. Supports env_chk.
            keys can be anything property of a kpoint object (kpts, kpts_shift,
            kpts_weights, labels, comment, coord_type, num_kpts,
            tet_connections, tet_number, tet_weight)
        input_filename (str): Input filename (if not "KPOINTS")
        output_filename (str): Output filename (if not "KPOINTS")
    """

    optional_params = [
        "kpoints_update",
        "input_filename",
        "output_filename",
    ]

    def run_task(self, fw_spec):

        kpoints_name = self.get("input_filename", "KPOINTS")
        kpoints = Kpoints.from_file(kpoints_name)

        kpoints_update = env_chk(self.get("kpoints_update"), fw_spec)

        if kpoints_update:
            for key, value in kpoints_update.items():
                setattr(kpoints, key, value)

        kpoints.write_file(self.get("output_filename", "KPOINTS"))


@explicit_serialize
class ModifyPotcar(FiretaskBase):
    """
    Modify Potcar file.

    Required params:
        potcar_symbols (dict): overwrite potcar with symbol. Supports env_chk.

    Optional params:
        functional (dict): functional to use, e.g. PBE, PBE_52, LDA_US, PW91
        input_filename (str): Input filename (if not "INCAR")
        output_filename (str): Output filename (if not "INCAR")
    """

    required_params = ["potcar_symbols"]
    optional_params = ["functional", "input_filename", "output_filename"]

    def run_task(self, fw_spec):
        potcar_symbols = env_chk(self.get("potcar_symbols"), fw_spec)
        functional = self.get("functional", None)
        potcar_name = self.get("input_filename", "POTCAR")
        potcar = Potcar.from_file(potcar_name)

        # Replace PotcarSingles corresponding to elements
        # contained in potcar_symbols
        for n, psingle in enumerate(potcar):
            if psingle.element in potcar_symbols:
                potcar[n] = PotcarSingle.from_symbol_and_functional(
                    potcar_symbols[psingle.element], functional
                )

        potcar.write_file(self.get("output_filename", "POTCAR"))


@explicit_serialize
class UpdateScanRelaxBandgap(FiretaskBase):
    """
    Writes input files for a SCAN relaxation by constructing a new input set.
    The purpose of this Firetask is to allow the KSPACING and smearing parameters
    to be recalculated based on the bandgap from the PBE relaxation in the
    SCAN relaxation workflow. Assumes that output files from a previous
    (e.g., optimization) run can be accessed in current dir or prev_calc_dir.

    Optional params (dict):
        override_default_vasp_params: Dict of any keyword arguments supported
                                      by MPScanRelaxSet.
        potcar_spec (bool): Instead of writing the POTCAR, write a
            "POTCAR.spec". This is intended to allow testing of workflows
            without requiring pseudo-potentials to be installed on the system.

    """
    optional_params = ["override_default_vasp_params", "potcar_spec"]

    def run_task(self, fw_spec):

        kwargs = self.get("override_default_vasp_params")
        potcar_spec = self.get("potcar_spec", False)

        os.chdir(os.getcwd())
        vrun = Vasprun("vasprun.xml", parse_potcar_file=False)
        bandgap = vrun.get_band_structure().get_band_gap()["energy"]
        structure = vrun.final_structure
        vis = MPScanRelaxSet(structure, bandgap=bandgap, **kwargs)
        vis.write_input(".", potcar_spec=potcar_spec)


@explicit_serialize
class WriteVaspStaticFromPrev(FiretaskBase):
    """
    Writes input files for a static run. Assumes that output files from a
    previous (e.g., optimization) run can be accessed in current dir or
    prev_calc_dir. Also allows lepsilon (dielectric constant) calcs.

    Optional params:
        potcar_spec (bool): Instead of writing the POTCAR, write a
            "POTCAR.spec". This is intended to allow testing of workflows
            without requiring pseudo-potentials to be installed on the system.
        (documentation for all other optional params can be found in
        MPStaticSet)

    """

    optional_params = [
        "prev_calc_dir",
        "reciprocal_density",
        "small_gap_multiply",
        "standardize",
        "sym_prec",
        "international_monoclinic",
        "lepsilon",
        "other_params",
        "potcar_spec",
    ]

    def run_task(self, fw_spec):
        lepsilon = self.get("lepsilon")

        # more k-points for dielectric calc.
        default_reciprocal_density = 200 if lepsilon else 100
        other_params = self.get("other_params", {})
        user_incar_settings = other_params.get("user_incar_settings", {})

        # for lepsilon runs, set EDIFF to 1E-5 unless user says otherwise
        if (
            lepsilon
            and "EDIFF" not in user_incar_settings
            and "EDIFF_PER_ATOM" not in user_incar_settings
        ):
            if "user_incar_settings" not in other_params:
                other_params["user_incar_settings"] = {}
            other_params["user_incar_settings"]["EDIFF"] = 1e-5

        vis = MPStaticSet.from_prev_calc(
            prev_calc_dir=self.get("prev_calc_dir", "."),
            reciprocal_density=self.get(
                "reciprocal_density", default_reciprocal_density
            ),
            small_gap_multiply=self.get("small_gap_multiply", None),
            standardize=self.get("standardize", False),
            sym_prec=self.get("sym_prec", 0.1),
            international_monoclinic=self.get(
                "international_monoclinic", True
            ),
            lepsilon=lepsilon,
            **other_params
        )

        potcar_spec = self.get("potcar_spec", False)
        vis.write_input(".", potcar_spec=potcar_spec)


@explicit_serialize
class WriteVaspHSEBSFromPrev(FiretaskBase):
    """
    Writes input files for HSE band structure run. Assumes that output files
    from a previous job can be accessed. Since HSE always re-optimizes the
    charge density (no nSCF mode), the previous job is used to get the location
    of VBM/CBM for mode="gap" (otherwise just used to get the structure /
    starting charge density).

    Optional params:
        potcar_spec (bool): Instead of writing the POTCAR, write a
            "POTCAR.spec". This is intended to allow testing of workflows
            without requiring pseudo-potentials to be installed on the system.
        (documentation for all other optional params can be found in
        MPHSEBSSet)
    """

    required_params = []
    optional_params = [
        "prev_calc_dir",
        "mode",
        "reciprocal_density",
        "kpoints_line_density",
        "potcar_spec",
    ]

    def run_task(self, fw_spec):
        vis = MPHSEBSSet.from_prev_calc(
            self.get("prev_calc_dir", "."),
            mode=self.get("mode", "uniform"),
            reciprocal_density=self.get("reciprocal_density", 50),
            kpoints_line_density=self.get("kpoints_line_density", 10),
            copy_chgcar=False,
        )
        potcar_spec = self.get("potcar_spec", False)
        vis.write_input(".", potcar_spec=potcar_spec)


@explicit_serialize
class WriteVaspNSCFFromPrev(FiretaskBase):
    """
    Writes input files for an NSCF static run. Assumes that output files from an
    scf job can be accessed. There are many options, e.g. uniform mode,
    line mode, adding the optical properties, etc.

    Optional params:
        potcar_spec (bool): Instead of writing the POTCAR, write a
            "POTCAR.spec". This is intended to allow testing of workflows
            without requiring pseudo-potentials to be installed on the system.
        (documentation for all optional params can be found in
        NonSCFVaspInputSet)
    """

    required_params = []
    optional_params = [
        "prev_calc_dir",
        "copy_chgcar",
        "nbands_factor",
        "reciprocal_density",
        "kpoints_line_density",
        "small_gap_multiply",
        "standardize",
        "sym_prec",
        "international_monoclinic",
        "mode",
        "nedos",
        "optics",
        "other_params",
        "potcar_spec",
    ]

    def run_task(self, fw_spec):
        vis = MPNonSCFSet.from_prev_calc(
            prev_calc_dir=self.get("prev_calc_dir", "."),
            copy_chgcar=self.get("copy_chgcar", False),
            nbands_factor=self.get("nbands_factor", 1.2),
            reciprocal_density=self.get("reciprocal_density", 100),
            kpoints_line_density=self.get("kpoints_line_density", 20),
            small_gap_multiply=self.get("small_gap_multiply", None),
            standardize=self.get("standardize", False),
            sym_prec=self.get("sym_prec", 0.1),
            international_monoclinic=self.get(
                "international_monoclinic", True
            ),
            mode=self.get("mode", "uniform"),
            nedos=self.get("nedos", 2001),
            optics=self.get("optics", False),
            **self.get("other_params", {})
        )
        potcar_spec = self.get("potcar_spec", False)
        vis.write_input(".", potcar_spec=potcar_spec)


@explicit_serialize
class WriteVaspSOCFromPrev(FiretaskBase):
    """
    Writes input files for a spinorbit coupling calculation.

    Required params:
        magmom (list): magnetic moment values for each site in the structure.
        saxis (list): magnetic field direction

    Optional params:
        potcar_spec (bool): Instead of writing the POTCAR, write a
            "POTCAR.spec". This is intended to allow testing of workflows
            without requiring pseudo-potentials to be installed on the system.
        (documentation for all optional params can be found in MPSOCSet)

    """

    required_params = ["magmom", "saxis"]

    optional_params = [
        "prev_calc_dir",
        "copy_chgcar",
        "nbands_factor",
        "reciprocal_density",
        "small_gap_multiply",
        "standardize",
        "sym_prec",
        "international_monoclinic",
        "other_params",
        "potcar_spec",
    ]

    def run_task(self, fw_spec):
        # TODO: @albalu - can saxis have a default value e.g. [001] and be an
        #  optional parameter? -computron
        # TODO: @albalu - can magmom be auto-parsed from the previous calc?
        #  -computron

        vis = MPSOCSet.from_prev_calc(
            prev_calc_dir=self.get("prev_calc_dir", "."),
            magmom=self["magmom"],
            saxis=self["saxis"],
            copy_chgcar=self.get("copy_chgcar", False),
            nbands_factor=self.get("nbands_factor", 1.2),
            reciprocal_density=self.get("reciprocal_density", 100),
            small_gap_multiply=self.get("small_gap_multiply", None),
            standardize=self.get("standardize", False),
            sym_prec=self.get("sym_prec", 0.1),
            international_monoclinic=self.get(
                "international_monoclinic", True
            ),
            **self.get("other_params", {})
        )
        potcar_spec = self.get("potcar_spec", False)
        vis.write_input(".", potcar_spec=potcar_spec)


@explicit_serialize
class WriteVaspNMRFromPrev(FiretaskBase):
    """
    Writes input files for a NMR calculation

    Optional params::
        prev_calc_dir: path to previous calculation, else current directory
        mode (str): the NMR calculation type: cs or efg, default is cs
        isotopes (list): list of isotopes to include, default is to include the
                         lowest mass quadrupolar isotope for all applicable
                         elements
        reciprocal_density (int): the reciprocal density for the kpoint mesh,
            defaults to 100
        other_params (dict) : any other params passed to MPNMRSet as a dict.
        potcar_spec (bool): Instead of writing the POTCAR, write a
            "POTCAR.spec". This is intended to allow testing of workflows
            without requiring pseudo-potentials to be installed on the system.
    """

    optional_params = [
        "prev_calc_dir",
        "mode",
        "isotopes",
        "reciprocal_density",
        "other_params",
        "potcar_spec",
    ]

    def run_task(self, fw_spec):
        vis = MPNMRSet.from_prev_calc(
            prev_calc_dir=self.get("prev_calc_dir", "."),
            mode=self.get("mode", "cs"),
            isotopes=self.get("isotopes", None),
            reciprocal_density=self.get("reciprocal_density", 100),
            **self.get("other_params", {})
        )
        potcar_spec = self.get("potcar_spec", False)
        vis.write_input(".", potcar_spec=potcar_spec)


@explicit_serialize
class WriteTransmutedStructureIOSet(FiretaskBase):
    """
    Apply the provided transformations to the input structure and write the
    input set for that structure. Reads structure from POSCAR if no structure
    provided. Note that if a transformation yields many structures from one,
    only the last structure in the list is used.

    Required params:
        structure (Structure): input structure
        transformations (list): list of names of transformation classes as
            defined in the modules in pymatgen.transformations
        vasp_input_set (VaspInputSet): VASP input set.

    Optional params:
        transformation_params (list): list of dicts where each dict specifies
            the input parameters to instantiate the transformation class in the
            transformations list.
        override_default_vasp_params (dict): additional user input settings.
        prev_calc_dir: path to previous calculation if using structure from
            another calculation.
        potcar_spec (bool): Instead of writing the POTCAR, write a
            "POTCAR.spec". This is intended to allow testing of workflows
            without requiring pseudo-potentials to be installed on the system.
    """

    required_params = ["structure", "transformations", "vasp_input_set"]
    optional_params = [
        "prev_calc_dir",
        "transformation_params",
        "override_default_vasp_params",
        "potcar_spec",
    ]

    def run_task(self, fw_spec):

        transformations = []
        transformation_params = self.get(
            "transformation_params",
            [{} for _ in range(len(self["transformations"]))],
        )
        for t in self["transformations"]:
            found = False
            t_cls = None
            for m in [
                "advanced_transformations",
                "defect_transformations",
                "site_transformations",
                "standard_transformations",
            ]:
                mod = import_module("pymatgen.transformations.{}".format(m))

                try:
                    t_cls = getattr(mod, t)
                    found = True
                    continue
                except AttributeError:
                    pass

            if not found:
                raise ValueError("Could not find transformation: {}".format(t))

            t_obj = t_cls(**transformation_params.pop(0))
            transformations.append(t_obj)

        # TODO: @matk86 - should prev_calc_dir use CONTCAR instead of POSCAR?
        #  Note that if current dir, maybe POSCAR is indeed best ... -computron
        structure = (
            self["structure"]
            if not self.get("prev_calc_dir", None)
            else Poscar.from_file(
                os.path.join(self["prev_calc_dir"], "POSCAR")
            ).structure
        )
        ts = TransformedStructure(structure)
        transmuter = StandardTransmuter([ts], transformations)
        final_structure = transmuter.transformed_structures[
            -1
        ].final_structure.copy()
        vis_orig = self["vasp_input_set"]
        vis_dict = vis_orig.as_dict()
        vis_dict["structure"] = final_structure.as_dict()
        vis_dict.update(self.get("override_default_vasp_params", {}) or {})
        vis = vis_orig.__class__.from_dict(vis_dict)

        potcar_spec = self.get("potcar_spec", False)
        vis.write_input(".", potcar_spec=potcar_spec)

        dumpfn(transmuter.transformed_structures[-1], "transformations.json")


@explicit_serialize
class WriteNormalmodeDisplacedPoscar(FiretaskBase):
    """
    Displace the structure from the previous calculation along the provided
    normal mode by the given amount and write the corresponding Poscar file.
    The fw_spec must contain a "normalmodes" key with "eigenvecs" sub-key that
    is likely produced by a previous calc.

    Required params:
        mode (int): normal mode index
        displacement (float): displacement along the normal mode in Angstroms
    """

    required_params = ["mode", "displacement"]

    def run_task(self, fw_spec):
        mode = self["mode"]
        disp = self["displacement"]
        structure = Structure.from_file("POSCAR")
        nm_eigenvecs = np.array(fw_spec["normalmodes"]["eigenvecs"])
        nm_norms = np.linalg.norm(nm_eigenvecs, axis=2)

        # displace the sites along the given normal mode: displacement vector
        # for each site = normalized eigen vector * amount of displacement
        nm_displacement = (
            nm_eigenvecs[mode, :, :] * disp / nm_norms[mode, :, np.newaxis]
        )
        for i, vec in enumerate(nm_displacement):
            structure.translate_sites(i, vec, frac_coords=False)

        # write the modified structure to poscar
        structure.to(fmt="poscar", filename="POSCAR")
