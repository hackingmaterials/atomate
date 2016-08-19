# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

from numpy.linalg import norm
from pymatgen.transformations.site_transformations import TranslateSitesTransformation

"""
This module defines tasks for writing vasp input sets for various types of vasp calculations
"""

import os
from six.moves import range

from fireworks import FireTaskBase, explicit_serialize
from fireworks.utilities.dict_mods import apply_mod

from pymatgen.alchemy.materials import TransformedStructure
from pymatgen.alchemy.transmuters import StandardTransmuter
from pymatgen.io.vasp import Incar, Poscar, Vasprun
from pymatgen.io.vasp.sets import MPStaticSet, MPNonSCFSet, MPSOCSet, MPHSEBSSet

from matmethods.utils.utils import env_chk

__author__ = 'Anubhav Jain, Shyue Ping Ong, Kiran Mathew'
__email__ = 'ajain@lbl.gov'


def load_class(mod, name):
    mod = __import__(mod, globals(), locals(), [name], 0)
    return getattr(mod, name)


@explicit_serialize
class WriteVaspFromIOSet(FireTaskBase):
    """
    Create VASP input files using implementations of pymatgen's
    AbstractVaspInputSet. An input set can be provided as an object or as a
    String/parameter combo.

    Required params:
        structure (Structure): structure
        vasp_input_set (AbstractVaspInputSet or str): Either a VaspInputSet
            object or a string name for the VASP input set (e.g.,
            "MPRelaxSet").

    Optional params:
        vasp_input_params (dict): When using a string name for VASP input set,
            use this as a dict to specify kwargs for instantiating the input
            set parameters. For example, if you want to change the
            user_incar_settings, you should provide:
            {"user_incar_settings": ...}.
            This setting is ignored if you provide the full object
            representation of a VaspInputSet rather than a String.
    """

    required_params = ["structure", "vasp_input_set"]
    optional_params = ["vasp_input_params"]

    def run_task(self, fw_spec):
        # if a full VaspInputSet object was provided
        if hasattr(self['vasp_input_set'], 'write_input'):
            vis = self['vasp_input_set']

        # if VaspInputSet String + parameters was provided
        else:
            vis_cls = load_class("pymatgen.io.vasp.sets",
                                 self["vasp_input_set"])
            vis = vis_cls(self["structure"],
                          **self.get("vasp_input_params", {}))

        vis.write_input(".")


@explicit_serialize
class WriteVaspFromPMGObjects(FireTaskBase):
    """
    Write VASP files using pymatgen objects.

    Required params:
        (none)

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
class ModifyIncar(FireTaskBase):
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

    optional_params = ["incar_update", "incar_multiply", "incar_dictmod",
                       "input_filename", "output_filename"]

    def run_task(self, fw_spec):

        # load INCAR
        incar_name = self.get("input_filename", "INCAR")
        incar = Incar.from_file(incar_name)

        # process FireWork env values via env_chk
        incar_update = env_chk(self.get('incar_update'), fw_spec)
        incar_multiply = env_chk(self.get('incar_multiply'), fw_spec)
        incar_dictmod = env_chk(self.get('incar_dictmod'), fw_spec)

        if incar_update:
            incar.update(incar_update)

        if incar_multiply:
            for k in incar_multiply:
                incar[k] = incar[k] * incar_multiply[k]

        if incar_dictmod:
            apply_mod(incar_dictmod, incar)

        # write INCAR
        incar.write_file(self.get("output_filename", "INCAR"))


@explicit_serialize
class WriteVaspStaticFromPrev(FireTaskBase):
    """
    Writes input files for a static run. Assumes that output files from a
    relaxation job can be accessed. Also allows lepsilon calcs.

    Required params:
        (none)

    Optional params:
        (documentation for all optional params can be found in
        MPStaticSet)
    """

    required_params = ["prev_calc_dir"]
    optional_params = ["reciprocal_density", "small_gap_multiply",
                       "standardize", "sym_prec", "international_monoclinic",
                       "lepsilon", "other_params"]

    def run_task(self, fw_spec):
        lepsilon = self.get("lepsilon")

        default_reciprocal_density = 100 if not lepsilon else 200
        other_params = self.get("other_params", {})

        # for lepsilon runs, set EDIFF to 1E-5 unless user says otherwise
        user_incar_settings = self.get("other_params", {}).get("user_incar_settings", {})

        if lepsilon and "EDIFF" not in user_incar_settings and "EDIFF_PER_ATOM" not in user_incar_settings:
            if "user_incar_settings" not in other_params:
                other_params["user_incar_settings"] = {}
            other_params["user_incar_settings"]["EDIFF"] = 1E-5

        vis = MPStaticSet.from_prev_calc(prev_calc_dir=self["prev_calc_dir"],
                                         reciprocal_density=default_reciprocal_density,
                                         small_gap_multiply=self.get("small_gap_multiply", None),
                                         standardize=self.get("standardize", False),
                                         sym_prec=self.get("sym_prec", 0.1),
                                         international_monoclinic=self.get("international_monoclinic", True),
                                         lepsilon=lepsilon, **other_params)
        vis.write_input(".")


@explicit_serialize
class WriteVaspHSEBSFromPrev(FireTaskBase):
    """
    Writes input files for HSE Gap run. Assumes that output files from a
    an NSCF job (for getting VBM/CBM) can be accessed.

    Required params:
        prev_calc_dir

    Optional params:
        (documentation for all optional params can be found in
        MPHSEBSSet)
    """

    required_params = ["prev_calc_dir"]
    optional_params = ["mode", "reciprocal_density"]

    def run_task(self, fw_spec):
        vis = MPHSEBSSet.from_prev_calc(self["prev_calc_dir"], mode=self.get("mode", "Uniform"),
                                        reciprocal_density=self.get("reciprocal_density", 50),
                                        copy_chgcar=False)
        vis.write_input(".")


@explicit_serialize
class WriteVaspNSCFFromPrev(FireTaskBase):
    """
    Writes input files for an NSCF static run. Assumes that output files from an
    scf job can be accessed. There are many options, e.g. uniform mode,
    line mode, adding the optical properties, etc.

    Required params:
        (none)

    Optional params:
        (documentation for all optional params can be found in
        NonSCFVaspInputSet)
    """

    required_params = ["prev_calc_dir"]
    optional_params = ["copy_chgcar", "nbands_factor", "reciprocal_density",
                       "kpoints_line_density", "small_gap_multiply",
                       "standardize", "sym_prec", "international_monoclinic",
                       "mode", "nedos", "optics", "other_params"]

    def run_task(self, fw_spec):
        vis = MPNonSCFSet.from_prev_calc(
            prev_calc_dir=self["prev_calc_dir"],
            copy_chgcar=self.get("copy_chgcar", False),
            nbands_factor=self.get("nbands_factor", 1.2),
            reciprocal_density=self.get("reciprocal_density", 100),
            kpoints_line_density=self.get("kpoints_line_density", 20),
            small_gap_multiply=self.get("small_gap_multiply", None),
            standardize=self.get("standardize", False),
            sym_prec=self.get("sym_prec", 0.1),
            international_monoclinic=self.get("international_monoclinic", True),
            mode=self.get("mode", "uniform"),
            nedos=self.get("nedos", 601),
            optics=self.get("optics", False),
            **self.get("other_params", {}))
        vis.write_input(".")


@explicit_serialize
class WriteVaspSOCFromPrev(FireTaskBase):
    """
    Writes input files for a spinorbit coupling calculation.

    Required params:
        prev_calc_dir: path to previous calculation
        magmom (list): magnetic moment values for each site in the structure.
        saxis (list): magnetic field direction

    Optional params:
        (none)
    """
    required_params = ["prev_calc_dir", "magmom", "saxis"]

    optional_params = ["copy_chgcar", "nbands_factor", "reciprocal_density",
                       "small_gap_multiply", "standardize", "sym_prec",
                       "international_monoclinic", "other_params"]

    def run_task(self, fw_spec):
        vis = MPSOCSet.from_prev_calc(
            prev_calc_dir=self["prev_calc_dir"],
            magmom=self["magmom"],
            saxis=self["saxis"],
            copy_chgcar=self.get("copy_chgcar", False),
            nbands_factor=self.get("nbands_factor", 1.2),
            reciprocal_density=self.get("reciprocal_density", 100),
            small_gap_multiply=self.get("small_gap_multiply", None),
            standardize=self.get("standardize", False),
            sym_prec=self.get("sym_prec", 0.1),
            international_monoclinic=self.get("international_monoclinic", True),
            **self.get("other_params", {}))
        vis.write_input(".")


@explicit_serialize
class WriteTransmutedStructureIOSet(FireTaskBase):
    """
    Apply the provided transformations to the input structure and write the
    input set for that structure. Reads structure from POSCAR if no structure provided

    Required params:
        structure (Structure): input structure
        transformations (list): list of names of transformation classes as defined in
            the modules in pymatgen.transformations
        vasp_input_set (string): string name for the VASP input set (e.g.,
            "MPRelaxSet").

    Optional params:
        transformation_params (list): list of dicts where each dict specifies
            the input parameters to instantiate the transformation class in
            the transformations list.
        vasp_input_params (dict): When using a string name for VASP input set,
            use this as a dict to specify kwargs for instantiating the input
            set parameters. For example, if you want to change the
            user_incar_settings, you should provide: {"user_incar_settings": ...}.
            This setting is ignored if you provide the full object
            representation of a VaspInputSet rather than a String.
        prev_calc_dir: path to previous calculation if using structure 
            from another calculation
    """

    required_params = ["structure", "transformations", "vasp_input_set"]
    optional_params = ["prev_calc_dir", "transformation_params", "vasp_input_params"]

    def run_task(self, fw_spec):

        vis_cls = load_class("pymatgen.io.vasp.sets", self["vasp_input_set"])

        transformations = []
        transformation_params = self.get("transformation_params",
                                         [{} for i in range(len(self["transformations"]))])
        for t in self["transformations"]:
            for m in ["advanced_transformations", "defect_transformations",
                      "site_transformations", "standard_transformations"]:
                mod = __import__("pymatgen.transformations." + m, globals(), locals(), [t], -1)
                try:
                    t_cls = getattr(mod, t)
                except AttributeError:
                    continue
                t_obj = t_cls(**transformation_params.pop(0))
                transformations.append(t_obj)

        structure = self['structure'] if 'prev_calc_dir' not in self else \
                Poscar.from_file(os.path.join(self['prev_calc_dir'], 'POSCAR')).structure
        ts = TransformedStructure(structure)
        transmuter = StandardTransmuter([ts], transformations)
        vis = vis_cls(transmuter.transformed_structures[-1].final_structure, **self.get("vasp_input_params", {}))
        vis.write_input(".")


@explicit_serialize
class WriteNormalmodeDisplacementIOSet(FireTaskBase):
    """
    Displace the structure from the previous calculation along the provided normal mode by the
    given amount and write the corresponding vasp input set for dielectric constant calculation.

    Required params:
        mode (int): normal mode index
        displacement (float): displacement along the normal mode in Angstroms
        vasp_input_set (DictVaspInputSet): vasp input set.

    Optional params:
        vasp_input_params (dict): user vasp input settings
    """

    required_params = ["mode", "displacement", "vasp_input_set"]
    optional_params = ["vasp_input_params"]

    def run_task(self, fw_spec):
        vrun = Vasprun('vasprun.xml.gz')
        structure = vrun.final_structure.copy()
        normalmode_eigenvecs = vrun.normalmode_eigenvecs
        nmodes, natoms, _ = normalmode_eigenvecs.shape
        # normalize the eigen vectors
        for i in range(nmodes):
            for j in range(natoms):
                normalmode_eigenvecs[i, j, :] = normalmode_eigenvecs[i, j, :] / norm(normalmode_eigenvecs[i, j, :])

        # displace the sites along the given normal mode
        normalmode_displacement = normalmode_eigenvecs[self["mode"], :, :] * self["displacement"]
        transformation = TranslateSitesTransformation(range(len(structure)), normalmode_displacement,
                                                      vector_in_frac_coords=False)
        ts = TransformedStructure(structure)
        transmuter = StandardTransmuter([ts], [transformation])

        # write the static vasp input set corresponding to the transmuted structure to compute epsilon
        vis = self["vasp_input_set"].__class__(transmuter.transformed_structures[-1].final_structure,
                                               lepsilon=True, **self.get("vasp_input_params", {}))
        vis.write_input(".")
