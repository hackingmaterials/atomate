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
    def write_input_from_prevrun(kpoints_density=1000, prev_dir=None, standardization_symprec=0.1, preserve_magmom=True, preserve_old_incar=True, output_dir="."):

        # get old structure, including MAGMOM decoration if desired
        structure = get_structure_from_prev_run(prev_dir, preserve_magmom=preserve_magmom)

        # standardize the structure if desired
        if standardization_symprec:
            sym_finder = SpacegroupAnalyzer(structure, symprec=standardization_symprec)
            structure = sym_finder.get_primitive_standard_structure()

        vis = StaticVaspInputSet(kpoints_density)
        vis.write_input(structure, output_dir)



        """
        if preserve_old_incar:
            # TODO: parse old incar
            # TODO: override MAGMOM and also the STATIC_SETTINGS in this INCAR
            pass

        # TODO: use old CHGCAR!! ICHG=1, if you are not standardizing

        # write incar
        # TODO: add options to MPStaticVaspInputset constructor
        mpsvip = MPStaticVaspInputSet()
        # TODO: get KPOINTS and write it
        # TODO: get POSCAR and write it
        # TODO: get POTCAR and write it
        # TODO: get INCAR and write it -- but keep in mind preserving the old INCAR

        # TODO: add an option to preserve the old KPOINTS??


        # TODO: see if anything is missing
        """






class MPStaticVaspInputSet(DictVaspInputSet):
    """
    Implementation of VaspInputSet overriding MaterialsProjectVaspInputSet
    for static calculations that typically follow relaxation runs.
    It is recommended to use the static from_previous_run method to construct
    the input set to inherit most of the functions.

    Args:
        kpoints_density (int): kpoints density for the reciprocal cell of
            structure. Might need to increase the default value when
            calculating metallic materials.
        sym_prec (float): Tolerance for symmetry finding

    kwargs:
        hubbard_off (bool): Whether to turn off Hubbard U if it is specified in
            config_dict ("MP Static"). Defaults to False, i.e., follow settings
            in config_dict.
        user_incar_settings (dict): User INCAR settings. This allows a user
            to override INCAR settings, e.g., setting a different MAGMOM for
            various elements or species.
        constrain_total_magmom (bool): Whether to constrain the total magmom
            (NUPDOWN in INCAR) to be the sum of the expected MAGMOM for all
            species. Defaults to False.
        sort_structure (bool): Whether to sort the structure (using the
            default sort order of electronegativity) before generating input
            files. Defaults to True, the behavior you would want most of the
            time. This ensures that similar atomic species are grouped
            together.
        ediff_per_atom (bool): Whether the EDIFF is specified on a per atom
            basis.
    """

    def __init__(self, kpoints_density=90, sym_prec=0.1, **kwargs):
        super(MPStaticVaspInputSet, self).__init__(
            "MP Static",
            loadfn(os.path.join(MODULE_DIR, "MPVaspInputSet.yaml")),
            **kwargs)
        self.incar_settings.update(
            {"IBRION": -1, "ISMEAR": -5, "LAECHG": True, "LCHARG": True,
             "LORBIT": 11, "LVHAR": True, "LWAVE": False, "NSW": 0,
             "ICHARG": 0, "EDIFF": 0.000001, "ALGO": "Normal"})
        self.kpoints_settings.update({"kpoints_density": kpoints_density})
        self.sym_prec = sym_prec

    def get_kpoints(self, structure, primitive_standard=False):
        """
        Get a KPOINTS file using the fully automated grid method. Uses
        Gamma centered meshes for hexagonal cells and Monk grids otherwise.

        Args:
            structure (Structure/IStructure): structure to get kpoints
            primitive_standard (Bool): whether the input structure is
            a primitive standardized cell
        """
        if not primitive_standard:
            structure = self.get_poscar(structure).structure
        self.kpoints_settings['grid_density'] = \
            self.kpoints_settings["kpoints_density"] * \
            structure.lattice.reciprocal_lattice.volume * \
            structure.num_sites
        return super(MPStaticVaspInputSet, self).get_kpoints(structure)

    def get_poscar(self, structure):
        """
        Get a POSCAR file with a primitive standardized cell of
        the giving structure.

        Args:
            structure (Structure/IStructure): structure to get POSCAR
        """
        sym_finder = SpacegroupAnalyzer(structure, symprec=self.sym_prec)
        return Poscar(sym_finder.get_primitive_standard_structure(False))

    @staticmethod
    def get_structure(vasp_run, outcar=None, initial_structure=False,
                      additional_info=False, sym_prec=0.1):
        """
        Process structure for static calculations from previous run.

        Args:
            vasp_run (Vasprun): Vasprun that contains the final structure
                from previous run.
            outcar (Outcar): Outcar that contains the magnetization info from
                previous run.
            initial_structure (bool): Whether to return the structure from
                previous run. Default is False.
            additional_info (bool):
                Whether to return additional symmetry info related to the
                structure. If True, return a list of the refined structure (
                conventional cell), the conventional standard structure,
                the symmetry dataset and symmetry operations of the
                structure (see SpacegroupAnalyzer doc for details).
            sym_prec (float): Tolerance for symmetry finding

        Returns:
            Returns the magmom-decorated structure that can be passed to get
            Vasp input files, e.g. get_kpoints.
        """
        if vasp_run.is_spin:
            if outcar and outcar.magnetization:
                magmom = {"magmom": [i['tot'] for i in outcar.magnetization]}
            else:
                magmom = {
                    "magmom": vasp_run.as_dict()['input']['parameters']
                    ['MAGMOM']}
        else:
            magmom = None
        structure = vasp_run.final_structure
        if magmom:
            structure = structure.copy(site_properties=magmom)
        sym_finder = SpacegroupAnalyzer(structure, symprec=sym_prec)
        if initial_structure:
            return structure
        elif additional_info:
            info = [sym_finder.get_refined_structure(),
                    sym_finder.get_conventional_standard_structure(False),
                    sym_finder.get_symmetry_dataset(),
                    sym_finder.get_symmetry_operations()]
            return [sym_finder.get_primitive_standard_structure(False),
                    info]
        else:
            return sym_finder.get_primitive_standard_structure(False)

    @staticmethod
    def from_previous_vasp_run(previous_vasp_dir, output_dir='.',
                               user_incar_settings=None,
                               make_dir_if_not_present=True,
                               kpoints_density=90, sym_prec=0.1):
        """
        Generate a set of Vasp input files for static calculations from a
        directory of previous Vasp run.

        Args:
            previous_vasp_dir (str): Directory containing the outputs(
                vasprun.xml and OUTCAR) of previous vasp run.
            output_dir (str): Directory to write the VASP input files for
                the static calculations. Defaults to current directory.
            make_dir_if_not_present (bool): Set to True if you want the
                directory (and the whole path) to be created if it is not
                present.
            kpoints_density (int): kpoints density for the reciprocal cell
                of structure. Might need to increase the default value when
                calculating metallic materials.
            sym_prec (float): Tolerance for symmetry finding
        """
        # Read input and output from previous run
        try:
            vasp_run = Vasprun(os.path.join(previous_vasp_dir, "vasprun.xml"),
                               parse_dos=False, parse_eigen=None)
            outcar = Outcar(os.path.join(previous_vasp_dir, "OUTCAR"))
            previous_incar = vasp_run.incar
            previous_kpoints = vasp_run.kpoints
        except:
            traceback.print_exc()
            raise RuntimeError("Can't get valid results from previous run. prev dir: {}".format(previous_vasp_dir))

        mpsvip = MPStaticVaspInputSet(kpoints_density=kpoints_density,
                                      sym_prec=sym_prec)
        structure = mpsvip.get_structure(vasp_run, outcar)

        mpsvip.write_input(structure, output_dir, make_dir_if_not_present)
        new_incar = mpsvip.get_incar(structure)

        # Use previous run INCAR and override necessary parameters
        previous_incar.update({"IBRION": -1, "ISMEAR": -5, "LAECHG": True,
                               "LCHARG": True, "LORBIT": 11, "LVHAR": True,
                               "LWAVE": False, "NSW": 0, "ICHARG": 0,
                               "ALGO": "Normal"})

        for incar_key in ["MAGMOM", "NUPDOWN"]:
            if new_incar.get(incar_key, None):
                previous_incar.update({incar_key: new_incar[incar_key]})
            else:
                previous_incar.pop(incar_key, None)

        # use new LDAUU when possible b/c the Poscar might have changed
        # representation
        if previous_incar.get('LDAU'):
            u = previous_incar.get('LDAUU', [])
            j = previous_incar.get('LDAUJ', [])
            if sum([u[x] - j[x] for x, y in enumerate(u)]) > 0:
                for tag in ('LDAUU', 'LDAUL', 'LDAUJ'):
                    previous_incar.update({tag: new_incar[tag]})
            # ensure to have LMAXMIX for GGA+U static run
            if "LMAXMIX" not in previous_incar:
                previous_incar.update({"LMAXMIX": new_incar["LMAXMIX"]})

        # Compare ediff between previous and staticinputset values,
        # choose the tighter ediff
        previous_incar.update({"EDIFF": min(previous_incar.get("EDIFF", 1),
                                            new_incar["EDIFF"])})

        # add user settings
        if user_incar_settings:
            previous_incar.update(user_incar_settings)
        previous_incar.write_file(os.path.join(output_dir, "INCAR"))

        # Perform checking on INCAR parameters
        if any([previous_incar.get("NSW", 0) != 0,
                previous_incar["IBRION"] != -1,
                previous_incar["LCHARG"] is not True,
               any([sum(previous_incar["LDAUU"]) <= 0,
                    previous_incar["LMAXMIX"] < 4])
               if previous_incar.get("LDAU") else False]):
            raise ValueError("Incompatible INCAR parameters!")

        # Prefer to use k-point scheme from previous run
        new_kpoints = mpsvip.get_kpoints(structure)
        if previous_kpoints.style != new_kpoints.style:
            if previous_kpoints.style == Kpoints.supported_modes.Monkhorst and \
                    SpacegroupAnalyzer(structure, 0.1).get_lattice_type() != \
                    "hexagonal":
                k_div = (kp + 1 if kp % 2 == 1 else kp
                         for kp in new_kpoints.kpts[0])
                Kpoints.monkhorst_automatic(k_div). \
                    write_file(os.path.join(output_dir, "KPOINTS"))
            else:
                Kpoints.gamma_automatic(new_kpoints.kpts[0]). \
                    write_file(os.path.join(output_dir, "KPOINTS"))
        else:
            new_kpoints.write_file(os.path.join(output_dir, "KPOINTS"))

