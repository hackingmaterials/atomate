# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import math
import os

from monty.os.path import zpath
from monty.serialization import loadfn

from matmethods.utils.utils import get_logger
from pymatgen.analysis.structure_matcher import StructureMatcher

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
    DEFAULT_SETTINGS = {"IBRION": -1, "ISMEAR": -5, "LAECHG": True, "LCHARG": True, "LORBIT": 11,
                        "LVHAR": True, "LVTOT": True, "LWAVE": False, "NSW": 0, "ICHARG": 0,
                        "EDIFF": 0.000001, "ALGO": "Fast"}

    def __init__(self, config_dict_override=None, reciprocal_density=100, **kwargs):
        d = kwargs
        d["name"] = "static"
        d["config_dict"] = loadfn(os.path.join(MODULE_DIR, "MPVaspInputSet.yaml"))
        d["force_gamma"] = True
        if "grid_density" in d["config_dict"]["KPOINTS"]:
            del d["config_dict"]["KPOINTS"]["grid_density"]
        d["config_dict"]["KPOINTS"]["reciprocal_density"] = reciprocal_density
        d["config_dict"]["INCAR"].update(self.DEFAULT_SETTINGS)
        for k in ["INCAR", "KPOINTS", "POSCAR", "POTCAR"]:
            if config_dict_override and config_dict_override.get(k):
                if k in d["config_dict"]:
                    d["config_dict"][k].update(config_dict_override[k])
                else:
                    d["config_dict"][k] = config_dict_override[k]
        super(StaticVaspInputSet, self).__init__(**d)

    @staticmethod
    def write_input_from_prevrun(config_dict_override=None, reciprocal_density=100,
                                 small_gap_multiply=None, prev_dir=None,
                                 standardization_symprec=0.1, international_monoclinic=True,
                                 preserve_magmom=True, preserve_old_incar=False, output_dir="."):
        """
        Args:
            config_dict_override (dict): like {"INCAR": {"NPAR": 2}} to
                                override the default input set
            reciprocal_density (int): density of k-mesh by reciprocal
                                    volume (defaults to 100)
            small_gap_multiply ([float, float]) - if the gap is less than 1st index,
                                multiply the default reciprocal_density by the 2nd index
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
            sym_finder = SpacegroupAnalyzer(structure, symprec=standardization_symprec)
            new_structure = sym_finder.get_primitive_standard_structure(
                international_monoclinic=international_monoclinic)
            logger.info("Validating new cell...")
            # the primitive structure finding has had several bugs in the past
            # defend through validation
            vpa_old = structure.volume/structure.num_sites
            vpa_new = new_structure.volume/new_structure.num_sites
            if abs(vpa_old - vpa_new)/vpa_old > 0.02:
                raise ValueError("Standardizing cell failed! VPA old: {}, VPA new: {}".format(vpa_old, vpa_new))

            sm = StructureMatcher()
            if not sm.fit(structure, new_structure):
                raise ValueError("Standardizing cell failed! Old structure doesn't match new.")

            structure = new_structure
            logger.info("Finished cell standardization procedure.")

        # multiply the reciprocal density if needed:
        if small_gap_multiply:
            prev_dir = prev_dir or os.curdir
            vasprun = Vasprun(zpath(os.path.join(prev_dir, "vasprun.xml")), parse_dos=False)
            gap = vasprun.eigenvalue_band_properties[0]
            if gap <= small_gap_multiply[0]:
                reciprocal_density = reciprocal_density * small_gap_multiply[1]

        vis = StaticVaspInputSet(config_dict_override=config_dict_override,
                                 reciprocal_density=reciprocal_density)
        vis.write_input(structure, output_dir)
        if preserve_old_incar:
            write_preserved_incar(vis, structure, prev_dir, config_dict_override, output_dir)


class NonSCFVaspInputSet(DictVaspInputSet):
    DEFAULT_SETTINGS = {"IBRION": -1, "ISMEAR": 0, "SIGMA": 0.001,
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
        d["config_dict"]["INCAR"].update(self.DEFAULT_SETTINGS)
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
                                 reciprocal_density=None, small_gap_multiply=None,
                                 prev_dir=None, mode="uniform", magmom_cutoff=0.1,
                                 nbands_factor=1.2, preserve_magmom=True,
                                 preserve_old_incar=False, output_dir="."):
        """
        Args:
            config_dict_override (dict): like {"INCAR": {"NPAR": 2}} to
                override the default input set
            reciprocal_density (int): density of k-mesh by reciprocal volume
                (defaults to 1000 in uniform, 20 in line mode)
            small_gap_multiply ([float, float]) - if the gap is less than 1st index,
                                multiply the default reciprocal_density by the 2nd index
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
        if reciprocal_density is None:
            reciprocal_density = 1000 if mode == "uniform" else 20

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

        # multiply the reciprocal density if needed for small gap compounds
        if small_gap_multiply:
            prev_dir = prev_dir or os.curdir
            vasprun = Vasprun(zpath(os.path.join(prev_dir, "vasprun.xml")), parse_dos=False)
            gap = vasprun.eigenvalue_band_properties[0]
            if gap <= small_gap_multiply[0]:
                reciprocal_density = reciprocal_density * small_gap_multiply[1]

        nscfvis = NonSCFVaspInputSet(config_dict_override=nscf_config_dict,
                                     reciprocal_density=reciprocal_density,
                                     mode=mode)
        nscfvis.write_input(structure, output_dir)
        if preserve_old_incar:
            write_preserved_incar(nscfvis, structure, config_dict_override,
                                  output_dir)


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
            incar = Incar.from_file(zpath(os.path.join(prev_dir, "INCAR")))
            if incar and incar.get("MAGMOM"):
                magmom = {"magmom": incar["MAGMOM"]}
            elif vasprun:
                magmom = {"magmom": vasprun.as_dict()['input']['parameters']['MAGMOM']}
            # dont trust the parsing mag stuff from outcar. last resort..
            elif outcar and outcar.magnetization:
                magmom = {"magmom": [i['tot'] for i in outcar.magnetization]}
            else:
                logger.warn("No MAGMOM found for the spin-polarized calculation !")
        else:
            magmom = None
        return structure.copy(site_properties=magmom)
    else:
        return Poscar.from_file(
            zpath(os.path.join(prev_dir, "CONTCAR"))).structure


def get_incar_from_prev_run(prev_dir, new_structure, default_settings=None,incar_dict_override=None):
    """
    Return incar object from the previous run with custom settings as well as
    structure dependent adjustments for parameters such as magmom and ldau.

    Args:
        new_structure (Structure): structure from the previous dir with or
            without modifications
        default_settings (dict): default settings
        prev_dir (str): path to the previous run directory
        incar_dict_override (dict): dictionary of Incar parameters to be
            overridden

    Returns:
        Incar object with default_settings overridden by previous incar and
        incar_dict_override, in that order
    """
    prev_incar = None
    try:
        prev_incar = Incar.from_file(zpath(os.path.join(prev_dir, "INCAR")))
        # the poscar is used only to get the ldau parameter mappings
        prev_poscar = Poscar.from_file(zpath(os.path.join(prev_dir, "POSCAR")))
    except:
        raise RuntimeError(
            "Can't get valid results from previous run. prev dir: {}".format(
                prev_dir))
    # apply the default parameter settings
    if default_settings:
        prev_incar.update(default_settings)
    # set structure dependent incar parameters
    # set MAGMOM
    # Note: site dependent parameters like MAGMOM need structure with the
    # corresponding site properties set.
    if prev_incar.get("MAGMOM", None):
        # new_incar obtained from get_incar will have the magmom settings(
        # based on the site properties('magmom', 'spin') of the structure)
        # set only if the default yaml input set has magmom set in it
        if hasattr(new_structure[0], "magmom"):
            mag = []
            for site in new_structure:
                mag.append(site.magmom)
            prev_incar["MAGMOM"] = mag
    # set LDAU parameters
    for incar_key in ["LDAUL", "LDAUU", "LDAUJ"]:
        if prev_incar.get(incar_key, None):
            new_poscar = Poscar(new_structure)
            set_params(incar_key, prev_incar, prev_poscar, new_poscar)
    # choose the tighter ediff
    tighter_ediff = min(prev_incar.get("EDIFF", 1), incar_dict_override["EDIFF"])
    # override the parameter settings
    if incar_dict_override:
        prev_incar.update(incar_dict_override)
    # set ediff
    prev_incar.update({"EDIFF": tighter_ediff})
    # sanity check
    prev_incar_dict = prev_incar.as_dict()
    check_list = []
    for k, v in default_settings.items():
        check_list.append(prev_incar_dict[k] != v)
    if any(check_list):
        raise ValueError("INCAR parameters not set properly!")
    return prev_incar


def get_param_mappings(param, incar, poscar):
    """
    Get the mapping for INCAR parameter 'param' for each atomic
    type in the poscar file from the values set in the provided incar.

    Args:
        param (str) : Atomic species dependent incar parameter
        incar (Incar): Incar object with the values for LDA+U parameters
            set.
        poscar (Poscar): Poscar object.

    Returns:
        mappings (dict): {param : {atomic_symbol: value}} format
    """
    mappings = {}
    if incar.get(param):
        mappings[param] = {}
        vals = incar[param]
        if isinstance(vals, list):
            for i, sym in enumerate(poscar.site_symbols):
                mappings[param][sym] = vals[i]
        elif isinstance(vals, dict):
            comp = poscar.structure.composition
            elements = sorted([el for el in comp.elements if comp[el] > 0],
                              key=lambda e: e.X)
            most_electroneg = elements[-1].symbol
            # if vals is in {"most_electroneg": {"symbol": value}} format,
            # the way they are set in the default inuputset yaml files for
            # LDAU
            if vals.get(most_electroneg,None):
                mappings[param] = vals[most_electroneg]
            else:
                mappings[param] = vals
        else:
            logger.error("Unknown format for specifying the values "
                         "for the parameter "
                         "{0}. Provided {1}".format(param, vals))
            raise ValueError
    return mappings


def set_params(param, incar, prev_poscar, new_poscar):
    """
    Set the parameter 'param' in the Incar file based on the atomic type
    info from the new_poscar and the parameter mappings obtained from
    the prev_poscar.

    Args:
        param (str): Atomic species dependent incar parameter
        incar (Incar): Incar object
        prev_poscar (Poscar): previous poscar
        new_poscar (Poscar): new poscar
    """
    mappings = get_param_mappings(param, incar, prev_poscar)
    if mappings:
        incar[param] = [mappings[param].get(sym,0) for sym in
                    new_poscar.site_symbols]


def write_preserved_incar(vis, structure, prev_dir, config_dict_override=None, output_dir="."):
    """
    Get the incar from the previous directory, update it based on the
    provided structure and input set and write it to the output
    directory.

    Args:
        vis (DictVaspInputSet): vasp input set
        structure (Structure): the structure that will be used to adjust
            incar parameters.
        prev_dir (str): directory from where the incar file settings will
            be taken.
        config_dict_override (dict): override settings.
        output_dir (str): directory where the incar file will be written.
    """
    new_incar = vis.get_incar(structure)
    incar_dict_override = config_dict_override.get("INCAR", {}) if config_dict_override else {}
    # choose tighter ediff
    incar_dict_override.update({"EDIFF": min(incar_dict_override.get("EDIFF", 1), new_incar.get("EDIFF", 1))})
    # incar from prev run
    incar = get_incar_from_prev_run(prev_dir, structure, vis.DEFAULT_SETTINGS,incar_dict_override=incar_dict_override)
    # set MAGMOM = 0 0 total_magnetization
    # use the previous values if not overridden
    # assumption previous calculation is non-collinear
    if not incar_dict_override.get("MAGMOM") and incar.get("LSORBIT") or incar.get("LNONCOLLINEAR"):
        val = []
        if incar.get("MAGMOM"):
            for m in incar["MAGMOM"]:
                # make sure the previous calc in non-collinear
                if not isinstance(m, list):
                    val.append([0, 0, m])
        incar["MAGMOM"] = val
    incar.write_file(os.path.join(output_dir, "INCAR"))
