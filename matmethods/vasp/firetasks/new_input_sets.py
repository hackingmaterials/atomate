import os

import math

from monty.os.path import zpath
from monty.serialization import loadfn

from pymatgen.io.vasp import Poscar, Vasprun, Outcar, Kpoints
from pymatgen.io.vasp.sets import DictVaspInputSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath

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
            vasprun = Vasprun(zpath(os.path.join(prev_dir, "vasprun.xml")), parse_dos=False, parse_eigen=False)
            outcar = Outcar(zpath(os.path.join(prev_dir, "OUTCAR")))
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
            return Poscar.from_file(zpath(os.path.join(prev_dir, "CONTCAR"))).structure


class StructureOptimizationVaspInputSet(DictVaspInputSet):

    def __init__(self, config_dict_override=None, reciprocal_density=50, force_gamma=True, **kwargs):
        d = kwargs
        d["name"] = "structure optimization"
        d["config_dict"] = loadfn(os.path.join(MODULE_DIR, "MPVaspInputSet.yaml"))
        if config_dict_override:
            d["config_dict"].update(config_dict_override)
        d["force_gamma"] = force_gamma
        if "grid_density" in d["config_dict"]["KPOINTS"]:
            del d["config_dict"]["KPOINTS"]["grid_density"]
        d["config_dict"]["KPOINTS"]["reciprocal_density"] = reciprocal_density

        super(StructureOptimizationVaspInputSet, self).__init__(**d)


class StaticVaspInputSet(DictVaspInputSet):

    STATIC_SETTINGS = {"IBRION": -1, "ISMEAR": -5, "LAECHG": True, "LCHARG": True,
             "LORBIT": 11, "LVHAR": True, "LVTOT": True, "LWAVE": False, "NSW": 0,
             "ICHARG": 0, "EDIFF": 0.000001, "ALGO": "Fast"}

    def __init__(self, config_dict_override=None, reciprocal_density=100, **kwargs):
        d = kwargs
        d["name"] = "static"
        d["config_dict"] = loadfn(os.path.join(MODULE_DIR, "MPVaspInputSet.yaml"))
        if config_dict_override:
            d["config_dict"].update(config_dict_override)
        d["force_gamma"] = True
        if "grid_density" in d["config_dict"]["KPOINTS"]:
            del d["config_dict"]["KPOINTS"]["grid_density"]
        d["config_dict"]["KPOINTS"]["reciprocal_density"] = reciprocal_density
        d["config_dict"]["INCAR"].update(self.STATIC_SETTINGS)
        super(StaticVaspInputSet, self).__init__(**d)

    @staticmethod
    def write_input_from_prevrun(config_dict_override=None, reciprocal_density=100, prev_dir=None, standardization_symprec=0.1, preserve_magmom=True, preserve_old_incar=False, output_dir="."):

        # get old structure, including MAGMOM decoration if desired
        structure = get_structure_from_prev_run(prev_dir, preserve_magmom=preserve_magmom)

        # standardize the structure if desired
        if standardization_symprec:
            sym_finder = SpacegroupAnalyzer(structure, symprec=standardization_symprec)
            structure = sym_finder.get_primitive_standard_structure()

        # TODO: re-use old CHGCAR (ICHG=1) if you are NOT standardizing (i.e., changing) the cell for faster performance. This is probably rare.

        vis = StaticVaspInputSet(config_dict_override=config_dict_override, reciprocal_density=reciprocal_density)
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
            # TODO: perform a final sanity check on the parameters(?)

        # TODO: add an option to preserve the old KPOINTS??


class NonSCFVaspInputSet(DictVaspInputSet):

    NSCF_SETTINGS = {"IBRION": -1, "ISMEAR": 0, "SIGMA": 0.001, "LCHARG": False,
             "LORBIT": 11, "LWAVE": False, "NSW": 0, "ISYM": 0, "ICHARG": 11}

    # TODO: kpoints density is not really correct
    def __init__(self, config_dict_override=None, mode="uniform", kpoints_density=None, sym_prec=0.1, **kwargs):
        if kpoints_density is None:
            kpoints_density = 1000 if mode == "uniform" else 20

        self.kpoints_density = kpoints_density  # used by the "get_kpoints()" method

        d = kwargs
        d["name"] = "non scf"
        d["config_dict"] = loadfn(os.path.join(MODULE_DIR, "MPVaspInputSet.yaml"))
        if config_dict_override:
            d["config_dict"].update(config_dict_override)
        d["config_dict"]["INCAR"].update(self.NSCF_SETTINGS)
        if mode == "uniform":
            d["config_dict"]["INCAR"].update({"NEDOS": 601})
        super(NonSCFVaspInputSet, self).__init__(**d)

        self.sym_prec = sym_prec
        self.mode = mode

    def get_kpoints(self, structure):
        """
        Get a KPOINTS file for NonSCF calculation. In "Line" mode, kpoints are
        generated along high symmetry lines. In "Uniform" mode, kpoints are
        Gamma-centered mesh grid. Kpoints are written explicitly in both cases.
        """

        if self.mode == "line":
            kpath = HighSymmKpath(structure)
            frac_k_points, k_points_labels = kpath.get_kpoints(line_density=self.kpoints_settings["kpoints_line_density"], coords_are_cartesian=False)
            return Kpoints(comment="Non SCF run along symmetry lines", style=Kpoints.supported_modes.Reciprocal,
                           num_kpts=len(frac_k_points), kpts=frac_k_points, labels=k_points_labels, kpts_weights=[1] * len(frac_k_points))
        else:
            num_kpoints = self.kpoints_settings["kpoints_density"] * structure.lattice.reciprocal_lattice.volume
            kpoints = Kpoints.automatic_density(structure, num_kpoints * structure.num_sites)
            mesh = kpoints.kpts[0]
            ir_kpts = SpacegroupAnalyzer(structure, symprec=self.sym_prec).get_ir_reciprocal_mesh(mesh)
            kpts = []
            weights = []
            for k in ir_kpts:
                kpts.append(k[0])
                weights.append(int(k[1]))
            return Kpoints(comment="Non SCF run on uniform grid", style=Kpoints.supported_modes.Reciprocal,
                           num_kpts=len(ir_kpts), kpts=kpts, kpts_weights=weights)

    @staticmethod
    def write_input_from_prevrun(mode="uniform", magmom_cutoff=0.2, nbands_factor=1.2, kpoints_density=None, prev_dir=None, preserve_magmom=True, preserve_old_incar=False, output_dir=".", config_dict_override=None):

        # TODO: the user_incar_settings are not used ... FIX THISSSSSS!!
        user_incar_settings = {}

        # get old structure, including MAGMOM decoration if desired
        structure = get_structure_from_prev_run(prev_dir, preserve_magmom=preserve_magmom)

        # crank up NBANDS by nbands_factor
        prev_dir = prev_dir or os.curdir
        vasprun = Vasprun(os.path.join(prev_dir, "vasprun.xml"), parse_dos=False, parse_eigen=False)
        user_incar_settings["NBANDS"] = int(math.ceil(vasprun.as_dict()["input"]["parameters"]["NBANDS"] * nbands_factor))

        # retain grid of old run
        for grid in ["NGX", "NGY", "NGZ"]:
            if vasprun.incar.get(grid):
                user_incar_settings[grid] = vasprun.incar.get(grid)

        if magmom_cutoff:
            # turn off ISPIN if previous calc did not have significant magnetic moments (>magmom_cutoff)
            if vasprun.is_spin:
                outcar = Outcar(zpath(os.path.join(prev_dir, "OUTCAR")))
                magmom_cutoff = [i['tot'] > magmom_cutoff for i in outcar.magnetization]
                ispin = 2 if any(magmom_cutoff) else 1
            else:
                ispin = 1
            user_incar_settings["ISPIN"] = ispin

        # TODO: add config_diect override
        nscfvis = NonSCFVaspInputSet(mode=mode, kpoints_density=kpoints_density)
        nscfvis.write_input(structure, output_dir)

        if preserve_old_incar:
            raise NotImplementedError("The option to preserve the old INCAR is not yet implemented!")
            # TODO: implement me!

        # TODO: perform final checks on parameters