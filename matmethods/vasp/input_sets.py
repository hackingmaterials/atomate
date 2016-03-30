# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import math
import os

from monty.os.path import zpath
from monty.serialization import loadfn

from matmethods.utils.utils import get_logger

from pymatgen.io.vasp import Incar, Poscar, Vasprun, Outcar, Kpoints
from pymatgen.io.vasp.sets import DictVaspInputSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath

__author__ = 'Anubhav Jain, Kiran Mathew'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov'

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

logger = get_logger(__name__)


class StructureOptimizationVaspInputSet(DictVaspInputSet):
    def __init__(self, config_dict_override=None, reciprocal_density=50,
                 force_gamma=True, **kwargs):
        d = kwargs
        d["name"] = "structure optimization"
        d["config_dict"] = loadfn(
            os.path.join(MODULE_DIR, "MPVaspInputSet.yaml"))

        for k in ["INCAR", "KPOINTS", "POSCAR", "POTCAR"]:
            if config_dict_override and config_dict_override.get(k):
                if k in d["config_dict"]:
                    d["config_dict"][k].update(config_dict_override[k])
                else:
                    d["config_dict"][k] = config_dict_override[k]
        d["force_gamma"] = force_gamma
        if "grid_density" in d["config_dict"]["KPOINTS"]:
            del d["config_dict"]["KPOINTS"]["grid_density"]
        d["config_dict"]["KPOINTS"]["reciprocal_density"] = reciprocal_density

        super(StructureOptimizationVaspInputSet, self).__init__(**d)


class StaticVaspInputSet(DictVaspInputSet):
    STATIC_SETTINGS = {"IBRION": -1, "ISMEAR": -5, "LAECHG": True,
                       "LCHARG": True,
                       "LORBIT": 11, "LVHAR": True, "LVTOT": True,
                       "LWAVE": False, "NSW": 0,
                       "ICHARG": 0, "EDIFF": 0.000001, "ALGO": "Fast"}

    def __init__(self, config_dict_override=None, reciprocal_density=100,
                 **kwargs):
        d = kwargs
        d["name"] = "static"
        d["config_dict"] = loadfn(
            os.path.join(MODULE_DIR, "MPVaspInputSet.yaml"))
        d["force_gamma"] = True
        if "grid_density" in d["config_dict"]["KPOINTS"]:
            del d["config_dict"]["KPOINTS"]["grid_density"]
        d["config_dict"]["KPOINTS"]["reciprocal_density"] = reciprocal_density
        d["config_dict"]["INCAR"].update(self.STATIC_SETTINGS)
        for k in ["INCAR", "KPOINTS", "POSCAR", "POTCAR"]:
            if config_dict_override and config_dict_override.get(k):
                if k in d["config_dict"]:
                    d["config_dict"][k].update(config_dict_override[k])
                else:
                    d["config_dict"][k] = config_dict_override[k]
        super(StaticVaspInputSet, self).__init__(**d)

    @staticmethod
    def write_input_from_prevrun(config_dict_override=None,
                                 reciprocal_density=100, prev_dir=None,
                                 standardization_symprec=0.1,
                                 international_monoclinic=True,
                                 preserve_magmom=True,
                                 preserve_old_incar=False, output_dir="."):
        """
        Args:
            config_dict_override (dict): like {"INCAR": {"NPAR": 2}} to
                                override the default input set
            reciprocal_density (int): density of k-mesh by reciprocal
                                    volume (defaults to 100)
            prev_dir(str): directory containing output files of the previous
                        relaxation run. Defaults to current dir.
            standardization_symprec (float): Symprec for standardization.
                        Set to None for no cell standardization. Defaults to 0.1.
            international_monoclinic (bool): Whether to use international
                    convention (vs Curtarolo) for monoclinic. Defaults True.
            preserve_magmom (bool): whether to preserve old MAGMOM. Defaults
                        to True
            preserve_old_incar (bool): whether to try to preserve most of
                            the older INCAR parameters instead of overriding
                            with Inputset values. Defaults False.
            output_dir (str): where to put the output files (defaults
                        current dir)
        """
        structure = get_structure_from_prev_run(prev_dir,
                                                preserve_magmom=preserve_magmom)
        # standardize the structure if desired
        if standardization_symprec:
            logger.info("Standardizing cell...")
            sym_finder = SpacegroupAnalyzer(structure,
                                            symprec=standardization_symprec)
            structure = sym_finder.get_primitive_standard_structure(
                international_monoclinic=international_monoclinic)
        vis = StaticVaspInputSet(config_dict_override=config_dict_override,
                                 reciprocal_density=reciprocal_density)
        vis.write_input(structure, output_dir)
        if preserve_old_incar:
            new_incar = vis.get_incar(structure)
            incar_dict_override = config_dict_override.get("INCAR", None) if \
                config_dict_override else None
            incar = get_incar_from_prev_run(new_incar, structure,
                                            StaticVaspInputSet.STATIC_SETTINGS,
                                            prev_dir,
                                            incar_dict_override=incar_dict_override)
            incar.write_file(os.path.join(output_dir, "INCAR"))


class NonSCFVaspInputSet(DictVaspInputSet):
    NSCF_SETTINGS = {"IBRION": -1, "ISMEAR": 0, "SIGMA": 0.001,
                     "LCHARG": False, "LORBIT": 11, "LWAVE": False,
                     "NSW": 0, "ISYM": 0, "ICHARG": 11}

    ALLOWED_MODES = ["line", "uniform"]

    def __init__(self, config_dict_override=None, mode="uniform",
                 reciprocal_density=None, sym_prec=0.1, **kwargs):

        if mode not in self.ALLOWED_MODES:
            raise ValueError(
                "{} is not an allowed 'mode'! Possible values are: {}".format(
                    mode, self.ALLOWED_MODES))

        if reciprocal_density is None:
            reciprocal_density = 1000 if mode == "uniform" else 20

        d = kwargs
        d["name"] = "non scf"
        d["config_dict"] = loadfn(
            os.path.join(MODULE_DIR, "MPVaspInputSet.yaml"))
        d["config_dict"]["INCAR"].update(self.NSCF_SETTINGS)
        if mode == "uniform":
            d["config_dict"]["INCAR"].update({"NEDOS": 601})
        for k in ["INCAR", "KPOINTS", "POSCAR", "POTCAR"]:
            if config_dict_override and config_dict_override.get(k):
                if k in d["config_dict"]:
                    d["config_dict"][k].update(config_dict_override[k])
                else:
                    d["config_dict"][k] = config_dict_override[k]
        super(NonSCFVaspInputSet, self).__init__(**d)
        # used by the "get_kpoints()" method
        self.reciprocal_density = reciprocal_density
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
            frac_k_points, k_points_labels = kpath.get_kpoints(
                line_density=self.reciprocal_density,
                coords_are_cartesian=False)
            return Kpoints(comment="Non SCF run along symmetry lines",
                           style=Kpoints.supported_modes.Reciprocal,
                           num_kpts=len(frac_k_points), kpts=frac_k_points,
                           labels=k_points_labels,
                           kpts_weights=[1] * len(frac_k_points))
        else:
            kpoints = Kpoints.automatic_density_by_vol(structure,
                                                       self.reciprocal_density,
                                                       force_gamma=True)
            mesh = kpoints.kpts[0]
            ir_kpts = SpacegroupAnalyzer(structure,
                                         symprec=self.sym_prec).get_ir_reciprocal_mesh(
                mesh)
            kpts = []
            weights = []
            for k in ir_kpts:
                kpts.append(k[0])
                weights.append(int(k[1]))
            return Kpoints(comment="Non SCF run on uniform grid",
                           style=Kpoints.supported_modes.Reciprocal,
                           num_kpts=len(ir_kpts), kpts=kpts,
                           kpts_weights=weights)

    @staticmethod
    def write_input_from_prevrun(config_dict_override=None,
                                 reciprocal_density=None, prev_dir=None,
                                 mode="uniform", magmom_cutoff=0.1,
                                 nbands_factor=1.2, preserve_magmom=True,
                                 preserve_old_incar=False, output_dir="."):
        """
        Args:
            config_dict_override (dict): like {"INCAR": {"NPAR": 2}} to
                override the default input set
            reciprocal_density (int): density of k-mesh by reciprocal volume
                (defaults to 1000 in uniform, 20 in line mode)
            prev_dir (str): directory containing output files of the
                previous relaxation run. Defaults to current dir.
            mode (str): either "uniform" (default) or "line"
            magmom_cutoff (float): if not None, ISPIN is turned off if all
                magnetic moments of previous run less than this
            nbands_factor (float): factor to increase NBANDS from previous run
            preserve_magmom (bool): whether to preserve old MAGMOM.
                Defaults to True
            preserve_old_incar (bool): whether to try to preserve most of
                the older INCAR parameters instead of overriding with Inputset
                values. Defaults to False.
            output_dir (str): where to put the output files (defaults
                current dir)
        """

        nscf_config_dict = {"INCAR": {}, "KPOINTS": {}}

        # get old structure, including MAGMOM decoration if desired
        structure = get_structure_from_prev_run(prev_dir,
                                                preserve_magmom=preserve_magmom)

        # crank up NBANDS by nbands_factor
        prev_dir = prev_dir or os.curdir
        vasprun = Vasprun(os.path.join(prev_dir, "vasprun.xml"),
                          parse_dos=False, parse_eigen=False)
        nscf_config_dict["INCAR"]["NBANDS"] = int(
            math.ceil(vasprun.as_dict()["input"]["parameters"][
                          "NBANDS"] * nbands_factor))
        # retain grid of old run
        for grid in ["NGX", "NGY", "NGZ"]:
            if vasprun.incar.get(grid):
                nscf_config_dict["INCAR"][grid] = vasprun.incar.get(grid)
        if magmom_cutoff:
            # turn off ISPIN if previous calc did not have significant
            # magnetic moments (>magmom_cutoff)
            if vasprun.is_spin:
                outcar = Outcar(zpath(os.path.join(prev_dir, "OUTCAR")))
                magmom_cutoff = [i['tot'] > magmom_cutoff for i in
                                 outcar.magnetization]
                ispin = 2 if any(magmom_cutoff) else 1
            else:
                ispin = 1
            nscf_config_dict["INCAR"]["ISPIN"] = ispin
        if config_dict_override:
            nscf_config_dict["INCAR"].update(config_dict_override["INCAR"])
            nscf_config_dict["KPOINTS"].update(config_dict_override["KPOINTS"])

        nscfvis = NonSCFVaspInputSet(config_dict_override=nscf_config_dict,
                                     reciprocal_density=reciprocal_density,
                                     mode=mode)
        nscfvis.write_input(structure, output_dir)
        if preserve_old_incar:
            new_incar = nscfvis.get_incar(structure)
            incar_dict_override = config_dict_override.get("INCAR", None) if \
                config_dict_override else None
            incar = get_incar_from_prev_run(new_incar, structure,
                                            NonSCFVaspInputSet.NSCF_SETTINGS,
                                            prev_dir,
                                            incar_dict_override=incar_dict_override)
            incar.write_file(os.path.join(output_dir, "INCAR"))


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
        vasprun = Vasprun(zpath(os.path.join(prev_dir, "vasprun.xml")),
                          parse_dos=False, parse_eigen=False)
        outcar = Outcar(zpath(os.path.join(prev_dir, "OUTCAR")))
        structure = vasprun.final_structure

        if vasprun.is_spin:
            if outcar and outcar.magnetization:
                magmom = {"magmom": [i['tot'] for i in outcar.magnetization]}
            else:
                magmom = {
                    "magmom": vasprun.as_dict()['input']['parameters'][
                        'MAGMOM']}
        else:
            magmom = None
        return structure.copy(site_properties=magmom)
    else:
        return Poscar.from_file(
            zpath(os.path.join(prev_dir, "CONTCAR"))).structure


def get_incar_from_prev_run(new_incar, new_structure, default_settings,
                            prev_dir,
                            incar_dict_override=None):
    """
    Return incar object from the previous run with custom settings as well as
    structure dependent adjustments for parameters such as magmom and ldau.
    Adapted from pymatgen.io.vasp.sets

    Args:
        new_incar (Incar): Incar from the new structure
        new_structure (Structure): structure from the previous dir with or
            without modifications
        default_settings (dict): default settings
        prev_dir (str): path to the previous run directory
        config_dict_override (dict): dictionary of Incar parameters to be
            overridden

    Returns:
        Incar object
    """
    prev_incar = None
    try:
        prev_incar = Incar.from_file(zpath(os.path.join(prev_dir, "INCAR")))
        prev_poscar = Poscar.from_file(zpath(os.path.join(prev_dir, "POSCAR")))
    except:
        raise RuntimeError(
            "Can't get valid results from previous run. prev dir: {}".format(
                prev_dir))
    # Override using the default parameter settings
    prev_incar.update(default_settings)
    for incar_key in ["MAGMOM", "NUPDOWN"]:
        if new_incar.get(incar_key, None):
            prev_incar.update({incar_key: new_incar[incar_key]})
        else:
            prev_incar.pop(incar_key, None)
    # set LDAU parameters
    if prev_incar.get('LDAU'):
        new_poscar = Poscar(new_structure)
        set_lda_params(prev_incar, prev_poscar, new_poscar)
        # ensure to have LMAXMIX for GGA+U static run
        if "LMAXMIX" not in prev_incar:
            prev_incar.update({"LMAXMIX": new_incar["LMAXMIX"]})
    # Choose the tighter ediff
    prev_incar.update(
        {"EDIFF": min(prev_incar.get("EDIFF", 1), new_incar["EDIFF"])})
    # override settings
    if incar_dict_override:
        prev_incar.update(incar_dict_override)
    # Sanity check
    prev_incar_dict = prev_incar.as_dict()
    check_list = []
    for k, v in default_settings.items():
        check_list.append(prev_incar_dict[k] != v)
    if any(check_list):
        raise ValueError("INCAR parameters not set properly!")
    return prev_incar


def get_lda_mappings(incar, poscar):
    """
    Get the lda+u parameters mapping for each atomic type in the poscar file
    from the values set in incar.

    Args:
        incar (Incar): Incar object with the values for LDA+U parameters
            set.
        poscar (Poscar): Poscar object.

    Returns:
        dict
    """
    lda_mappings = {}
    for lda_param in ("LDAUL", "LDAUU", "LDAUJ"):
        if incar.get(lda_param):
            vals = incar[lda_param]
            if isinstance(vals, list):
                lda_mappings[lda_param] = {}
                for i, sym in enumerate(poscar.site_symbols):
                    lda_mappings[lda_param][sym] = vals[i]
    return lda_mappings


def set_lda_params(incar, prev_poscar, new_poscar):
    """
    Set the LDA+U parameters  in the Incar file based on the atomic type
    info from the new_poscar and the ladau parameter mappings obtained from
    the prev_poscar.

    Args:
        incar (Incar): Incar object
        prev_poscar (Poscar): previous poscar
        new_poscar (Poscar): new poscar
    """
    lda_mappings = get_lda_mappings(incar, prev_poscar)
    for param, mappings in lda_mappings.items():
        incar[param] = [mappings[sym] for sym in new_poscar.site_symbols]
