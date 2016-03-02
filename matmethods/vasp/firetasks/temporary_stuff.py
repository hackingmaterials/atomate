import os

import math
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
            vasprun = Vasprun(os.path.join(prev_dir, "vasprun.xml"), parse_dos=False, parse_eigen=False)
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
            # TODO: perform a final sanity check on the parameters(?)

        # TODO: add an option to preserve the old KPOINTS??


class NonSCFVaspInputSet(DictVaspInputSet):

    NSCF_SETTINGS = {"IBRION": -1, "ISMEAR": 0, "SIGMA": 0.001, "LCHARG": False,
             "LORBIT": 11, "LWAVE": False, "NSW": 0, "ISYM": 0, "ICHARG": 11}

    # TODO: kpoints density is not really correct
    def __init__(self, mode="uniform", kpoints_density=None, sym_prec=0.1, **kwargs):
        if kpoints_density is None:
            kpoints_density = 1000 if mode == "uniform" else 20

        super(NonSCFVaspInputSet, self).__init__("MP Non-SCF",
            loadfn(os.path.join(MODULE_DIR, "MPVaspInputSet.yaml")), **kwargs)

        self.incar_settings.update(self.NSCF_SETTINGS)
        self.sym_prec = sym_prec
        self.mode = mode

        if mode == "uniform":
            self.incar_settings.update({"NEDOS": 601})
            self.kpoints_settings.update({"kpoints_density": kpoints_density})
        elif mode == "line":
            self.kpoints_settings.update({"kpoints_line_density": kpoints_density})
        else:
            raise ValueError("Supported modes for NonSCF runs are 'line' and 'uniform'!")

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
    def write_input_from_prevrun(mode="uniform", magmom_cutoff=0.2, nbands_factor=1.2, kpoints_density=None, prev_dir=None, preserve_magmom=True, preserve_old_incar=False, output_dir=".", user_incar_settings=None):

        user_incar_settings = user_incar_settings or {}

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
                outcar = Outcar(os.path.join(prev_dir, "OUTCAR"))
                magmom_cutoff = [i['tot'] > magmom_cutoff for i in outcar.magnetization]
                ispin = 2 if any(magmom_cutoff) else 1
            else:
                ispin = 1
            user_incar_settings["ISPIN"] = ispin

        nscfvis = NonSCFVaspInputSet(mode=mode, kpoints_density=kpoints_density, user_incar_settings=user_incar_settings)
        nscfvis.write_input(structure, output_dir)

        if preserve_old_incar:
            raise NotImplementedError("The option to preserve the old INCAR is not yet implemented!")
            # TODO: implement me!

        # TODO: perform final checks on parameters