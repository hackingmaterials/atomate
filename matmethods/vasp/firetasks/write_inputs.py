# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
This module defines tasks for writing vasp input sets for various types of
vasp calculations
"""

from fireworks import FireTaskBase, explicit_serialize
from fireworks.utilities.dict_mods import apply_mod
from pymatgen.io.vasp import Incar

from matmethods.utils.utils import env_chk
from matmethods.vasp.input_sets import StaticVaspInputSet, \
    NonSCFVaspInputSet

__author__ = 'Anubhav Jain, Shyue Ping Ong'
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
            "MPVaspInputSet").

    Optional params:
        vasp_input_params (dict): When using a string name for VASP input set,
            use this as a dict to specify kwargs for instantiating the input
            set parameters. For example, if you want to change the
            user_incar_settings, you should provide: {"user_incar_settings": ...}.
            This setting is ignored if you provide the full object
            representation of a VaspInputSet rather than a String.
    """

    required_params = ["structure", "vasp_input_set"]
    optional_params = ["vasp_input_params"]

    def run_task(self, fw_spec):
        structure = self['structure']

        # if a full VaspInputSet object was provided
        if hasattr(self['vasp_input_set'], 'write_input'):
            vis = self['vasp_input_set']

        # if VaspInputSet String + parameters was provided
        else:
            mod = __import__("pymatgen.io.vasp.sets", globals(), locals(),
                             [self["vasp_input_set"]], -1)
            vis = load_class("pymatgen.io.vasp.sets", self["vasp_input_set"])(
                **self.get("vasp_input_params", {}))

        vis.write_input(structure, ".")


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
    relaxation job can be accessed.

    Required params:
        (none)

    Optional params:
        (documentation for all optional params can be found in
        StaticVaspInputSet.write_input_from_prevrun)
    """

    optional_params = ["config_dict_override", "reciprocal_density",
                       "small_gap_multiply", "prev_dir", "standardization_symprec",
                       "international_monoclinic", "preserve_magmom",
                       "preserve_old_incar"]

    def run_task(self, fw_spec):
        StaticVaspInputSet.write_input_from_prevrun(
            config_dict_override=self.get("config_dict_override"),
            reciprocal_density=self.get("reciprocal_density", 100),
            small_gap_multiply=self.get("small_gap_multiply", None),
            prev_dir=self.get("prev_dir"),
            standardization_symprec=self.get("standardization_symprec", 0.1),
            international_monoclinic=self.get("international_monoclinic",
                                              True),
            preserve_magmom=self.get("preserve_magmom", True),
            preserve_old_incar=self.get("preserve_old_incar", False))


@explicit_serialize
class WriteVaspNSCFFromPrev(FireTaskBase):
    """
    Writes input files for a static run. Assumes that output files from an
    scf job can be accessed.

    Required params:
        (none)

    Optional params:
        (documentation for all optional params can be found in
        NonSCFVaspInputSet.write_input_from_prevrun)
    """

    optional_params = ["config_dict_override", "reciprocal_density",
                       "small_gap_multiply", "prev_dir", "mode", "magmom_cutoff",
                       "nbands_factor", "preserve_magmom",
                       "preserve_old_incar"]

    def run_task(self, fw_spec):
        NonSCFVaspInputSet.write_input_from_prevrun(
            config_dict_override=self.get("config_dict_override"),
            reciprocal_density=self.get("reciprocal_density"),
            small_gap_multiply=self.get("small_gap_multiply", None),
            prev_dir=self.get("prev_dir"),
            mode=self.get("mode", "uniform"),
            magmom_cutoff=self.get("magmom_cutoff", 0.1),
            nbands_factor=self.get("nbands_factor", 1.2),
            preserve_magmom=(self.get("preserve_magmom"), True),
            preserve_old_incar=self.get("preserve_old_incar", False))
