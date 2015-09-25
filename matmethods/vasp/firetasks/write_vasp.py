from fireworks import FireTaskBase, explicit_serialize
from fireworks.utilities.dict_mods import apply_mod
from pymatgen.io.vasp import Incar

__author__ = 'Anubhav Jain <ajain@lbl.gov>, Shyue Ping Ong <ongsp@ucsd.edu>'

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

    def _load_class(mod, name):
        mod = __import__(mod, globals(), locals(), [name], 0)
        return getattr(mod, name)

    def run_task(self, fw_spec):
        structure = self['structure']

        # if a full VaspInputSet object was provided
        if hasattr(self['vasp_input_set'], 'write_input'):
            vis = self['vasp_input_set']

        # if VaspInputSet String + parameters was provided
        else:
            mod = __import__("pymatgen.io.vasp.sets", globals(), locals(),
                             [self["vasp_input_set"]], -1)
            vis = self._load_class("pymatgen.io.vasp.sets", self["vasp_input_set"])(
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
    Modify an INCAR file

    Required params:
        (none)

    Optional params:
        key_update (dict): overwrite Incar dict key
        key_multiply ({<str>:<float>}) - multiply Incar key by a constant factor
        key_dictmod ([{}]): use DictMod language to change Incar
        incar_filename (str): Incar filename (if not "INCAR")
    """

    optional_params = ["key_update", "key_multiply", "key_dictmod", "incar_filename"]

    def run_task(self, fw_spec):
        incar_name = self.get("incar_filename", "INCAR")

        incar = Incar.from_file(incar_name)

        if 'key_update' in self:
            incar.update(self['key_update'])

        if 'key_multiply' in self:
            for k in self['key_multiply']:
                incar[k] = incar[k] * self['key_multiply'][k]

        if "key_dictmod" in self:
            apply_mod(self['key_dictmod'], incar)

        incar.write_file(incar_name)