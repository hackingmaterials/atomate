from fireworks import FireTaskBase, explicit_serialize
from fireworks.utilities.dict_mods import apply_mod
from matmethods.utils.utils import env_chk
from matmethods.vasp.firetasks.temporary_stuff import StaticVaspInputSet
from pymatgen.io.vasp import Incar

__author__ = 'Anubhav Jain <ajain@lbl.gov>, Shyue Ping Ong <ongsp@ucsd.edu>'


def load_class(mod, name):
    mod = __import__(mod, globals(), locals(), [name], 0)
    return getattr(mod, name)


@explicit_serialize
class WriteVaspFromIOSet(FireTaskBase):
    """
    Create VASP input files using implementations of pymatgen's AbstractVaspInputSet.
    An input set can be provided as an object or as a String/parameter combo.

    Required params:
        structure (Structure): structure
        vasp_input_set (AbstractVaspInputSet or str): Either a VaspInputSet object or
            a string name for the VASP input set (e.g., "MPVaspInputSet").

    Optional params:
        vasp_input_params (dict): When using a string name for VASP input set, use
            this as a dict to specify kwargs for instantiating the input set
            parameters. For example, if you want to change the user_incar_settings, you
            should provide: {"user_incar_settings": ...}. This setting is ignored if you
            provide the full object representation of a VaspInputSet rather than a String.
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
        key_update (dict): overwrite Incar dict key. Supports env_chk.
        key_multiply ([{<str>:<float>}]) - multiply Incar key by a constant factor. Supports env_chk.
        key_dictmod ([{}]): use DictMod language to change Incar. Supports env_chk.
        input_filename (str): Input filename (if not "INCAR")
        output_filename (str): Output filename (if not "INCAR")
    """

    optional_params = ["key_update", "key_multiply", "key_dictmod", "input_filename", "output_filename"]

    def run_task(self, fw_spec):

        # load INCAR
        incar_name = self.get("input_filename", "INCAR")
        incar = Incar.from_file(incar_name)

        # process FireWork env values via env_chk
        key_update = env_chk(self.get('key_update'), fw_spec)
        key_multiply = env_chk(self.get('key_multiply'), fw_spec)
        key_dictmod = env_chk(self.get('key_dictmod'), fw_spec)

        if key_update:
            incar.update(key_update)

        if key_multiply:
            for k in key_multiply:
                incar[k] = incar[k] * key_multiply[k]

        if key_dictmod:
            apply_mod(key_dictmod, incar)

        # write INCAR
        incar.write_file(self.get("output_filename", "INCAR"))


@explicit_serialize
class WriteVaspStaticFromPrev(FireTaskBase):
    """
    Writes input files for a static run. Assumes that output files from a relaxation job can be accessed.

    Required params:
        (none)

    Optional params:
        prev_dir (str): directory containing output files of the previous relaxation run. Defaults to current dir.
        standardization_symprec (float): Symprec for standardization. Set to None for no cell standardization. Defaults to 0.1.
        preserve_magmom (bool): whether to preserve old MAGMOM. Defaults to True
        preserve_old_incar (bool): whether to try to preserve most of the older INCAR parameters instead of overriding with Inputset values. Defaults to False.
    """

    optional_params = ["standardization_symprec", "prev_dir", "preserve_magmom", "preserve_old_incar"]

    def run_task(self, fw_spec):
        StaticVaspInputSet.write_input_from_prevrun(prev_dir=self.get("prev_dir", None),
                                                    standardization_symprec=self.get("standardization_symprec"),
                                                    preserve_magmom=self.get("preserve_magmom", True),
                                                    preserve_old_incar=self.get("preserve_old_incar", False))


class WriteVaspUniformFromPrev(FireTaskBase):
    """
    Writes input files for a static run. Assumes that output files from an scf job can be accessed.

    Required params:
        (none)

    Optional params:
        prev_dir (str): directory containing output files of the previous relaxation run. Defaults to current dir.
        standardization_symprec (float): Symprec for standardization. Set to None for no cell standardization. Defaults to 0.1.
        preserve_magmom (bool): whether to preserve old MAGMOM. Defaults to True
        preserve_old_incar (bool): whether to try to preserve most of the older INCAR parameters instead of overriding with Inputset values. Defaults to False.
    """

    optional_params = ["standardization_symprec", "prev_dir", "preserve_magmom", "preserve_old_incar"]

    def run_task(self, fw_spec):
        StaticVaspInputSet.write_input_from_prevrun(prev_dir=self.get("prev_dir", None),
                                                    standardization_symprec=self.get("standardization_symprec"),
                                                    preserve_magmom=self.get("preserve_magmom", True),
                                                    preserve_old_incar=self.get("preserve_old_incar", False))