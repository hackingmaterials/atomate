import os
import traceback

from monty.serialization import loadfn

from pymatgen.io.vasp import Poscar, Vasprun, Outcar, Kpoints
from pymatgen.io.vasp.sets import DictVaspInputSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

def get_structure_from_prev_run(prev_dir, preserve_magmom=True):
        """
        Process structure for static calculations from previous run.

        Args:
            prev_dir (str): directory of the previous run
            preserve_magmom (bool): whether to preserve magmom of the old run.

        Returns:
            Returns the magmom-decorated structure.
        """

        prev_dir = prev_dir or os.curdir

        if preserve_magmom:
            vasprun = Vasprun(os.path.join(prev_dir, "vasprun.xml"))
            outcar = Outcar(os.path.join(prev_dir, "OUTCAR"))
            structure = vasprun.final_structure

            if vasprun.is_spin:
                if outcar and outcar.magnetization:
                    magmom = {"magmom": [i['tot'] for i in outcar.magnetization]}
                else:
                    magmom = {
                        "magmom": vasprun.as_dict()['input']['parameters']['MAGMOM']}
            else:
                magmom = None
            return structure.copy(site_properties=magmom)
        else:
            return Poscar.from_file(os.path.join(prev_dir, "CONTCAR")).structure


class StaticVaspInputSet(DictVaspInputSet):

    STATIC_SETTINGS = {"IBRION": -1, "ISMEAR": -5, "LAECHG": True, "LCHARG": True,
             "LORBIT": 11, "LVHAR": True, "LVTOT": True, "LWAVE": False, "NSW": 0,
             "ICHARG": 0, "EDIFF": 0.000001, "ALGO": "Fast"}

    # TODO: kpoints density is not really correct
    def __init__(self, kpoints_density=1000, **kwargs):
        super(StaticVaspInputSet, self).__init__("MP Static",
            loadfn(os.path.join(MODULE_DIR, "MPVaspInputSet.yaml")), **kwargs)

        self.incar_settings.update(self.STATIC_SETTINGS)
        self.kpoints_settings.update({"kpoints_density": kpoints_density})

    @staticmethod
    def write_input_from_prevrun(kpoints_density=1000, prev_dir=None, standardization_symprec=0.1, preserve_magmom=True, preserve_old_incar=False, output_dir=".", user_incar_settings=None):

        user_incar_settings = user_incar_settings or {}

        # get old structure, including MAGMOM decoration if desired
        structure = get_structure_from_prev_run(prev_dir, preserve_magmom=preserve_magmom)

        # standardize the structure if desired
        if standardization_symprec:
            sym_finder = SpacegroupAnalyzer(structure, symprec=standardization_symprec)
            structure = sym_finder.get_primitive_standard_structure()

        # TODO: re-use old CHGCAR (ICHG=1) if you are not standardizing (i.e., changing) the cell.

        vis = StaticVaspInputSet(kpoints_density, user_incar_settings=user_incar_settings)
        vis.write_input(structure, output_dir)

        if preserve_old_incar:
            raise NotImplementedError("The option to preserve the old INCAR is not yet implemented!")
            # TODO: parse old incar
            # TODO: override STATIC_SETTINGS in this INCAR
            # TODO: make sure MAGMOM aligns correctly with sites in newest INCAR
            # TODO: make sure LDAU aligns correctly with sites in newest INCAR
            # TODO: make sure to use the tighter EDIFF
            # TODO: write the new INCAR
            # TODO: check old code to see if anything needed is missing

        # TODO: add an option to preserve the old KPOINTS??


class NonSCFVaspInputSet(DictVaspInputSet):

    NSCF_SETTINGS = {"IBRION": -1, "ISMEAR": 0, "SIGMA": 0.001, "LCHARG": False,
             "LORBIT": 11, "LWAVE": False, "NSW": 0, "ISYM": 0, "ICHARG": 11}

    # TODO: kpoints density is not really correct
    def __init__(self, mode="uniform", kpoints_density=1000, **kwargs):
        if mode not in ["line", "uniform"]:
            raise ValueError("Supported modes for NonSCF runs are 'line' and 'uniform'!")

        super(NonSCFVaspInputSet, self).__init__("MP Non-SCF",
            loadfn(os.path.join(MODULE_DIR, "MPVaspInputSet.yaml")), **kwargs)

        self.incar_settings.update(self.NSCF_SETTINGS)

        self.kpoints_settings.update({"kpoints_density": kpoints_density})
        if mode == "uniform":
            self.incar_settings.update({"NEDOS": 601})

        if "NBANDS" not in kwargs:
            raise KeyError("For NonSCF runs, NBANDS value from SC runs is required!")

        # self.incar_settings.update(user_incar_settings)  TODO: is this needed?




    @staticmethod
    def write_input_from_prevrun(kpoints_density=1000, prev_dir=None, standardization_symprec=0.1, preserve_magmom=True, preserve_old_incar=False, output_dir=".", user_incar_settings=None):

        user_incar_settings = user_incar_settings or {}

        # get old structure, including MAGMOM decoration if desired
        structure = get_structure_from_prev_run(prev_dir, preserve_magmom=preserve_magmom)

        # standardize the structure if desired
        if standardization_symprec:
            sym_finder = SpacegroupAnalyzer(structure, symprec=standardization_symprec)
            structure = sym_finder.get_primitive_standard_structure()

        # TODO: re-use old CHGCAR (ICHG=1) if you are not standardizing (i.e., changing) the cell.

        vis = StaticVaspInputSet(kpoints_density, user_incar_settings=user_incar_settings)
        vis.write_input(structure, output_dir)

        if preserve_old_incar:
            raise NotImplementedError("The option to preserve the old INCAR is not yet implemented!")
            # TODO: parse old incar
            # TODO: override STATIC_SETTINGS in this INCAR
            # TODO: make sure MAGMOM aligns correctly with sites in newest INCAR
            # TODO: make sure LDAU aligns correctly with sites in newest INCAR
            # TODO: make sure to use the tighter EDIFF
            # TODO: write the new INCAR
            # TODO: check old code to see if anything needed is missing

        # TODO: add an option to preserve the old KPOINTS??