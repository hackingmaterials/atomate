from pymatgen.io.vasp.sets import DictVaspInputSet


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



class MPNonSCFVaspInputSet(MPStaticVaspInputSet):
    """
    Implementation of VaspInputSet overriding MaterialsProjectVaspInputSet
    for non self-consistent field (NonSCF) calculation that follows
    a static run to calculate bandstructure, density of states(DOS) and etc.
    It is recommended to use the NonSCF from_previous_run method to construct
    the input set to inherit most of the functions.

    Args:
        user_incar_settings (dict): A dict specify customized settings
            for INCAR. Must contain a NBANDS value, suggest to use
            1.2*(NBANDS from static run).
        mode: Line: Generate k-points along symmetry lines for
            bandstructure. Uniform: Generate uniform k-points
            grids for DOS.
        constrain_total_magmom (bool): Whether to constrain the total
            magmom (NUPDOWN in INCAR) to be the sum of the expected
            MAGMOM for all species. Defaults to False.
        kpoints_density (int): kpoints density for the reciprocal cell
            of structure. Might need to increase the default value when
            calculating metallic materials.
        kpoints_line_density (int): kpoints density to use in line-mode.
            Might need to increase the default value when calculating
            metallic materials.
        sort_structure (bool): Whether to sort structure. Defaults to
            False.
        sym_prec (float): Tolerance for symmetry finding
    """

    def __init__(self, user_incar_settings, mode="Line",
                 constrain_total_magmom=False, sort_structure=False,
                 kpoints_density=1000, sym_prec=0.1, kpoints_line_density=20):
        self.mode = mode
        self.sym_prec = sym_prec
        self.kpoints_line_density = kpoints_line_density
        if mode not in ["Line", "Uniform"]:
            raise ValueError("Supported modes for NonSCF runs are 'Line' and "
                             "'Uniform'!")
        DictVaspInputSet.__init__(self,
            "Materials Project Static",
            loadfn(os.path.join(MODULE_DIR, "MPVaspInputSet.yaml")),
            constrain_total_magmom=constrain_total_magmom,
            sort_structure=sort_structure)
        self.user_incar_settings = user_incar_settings
        self.incar_settings.update(
            {"IBRION": -1, "ISMEAR": 0, "SIGMA": 0.001, "LCHARG": False,
             "LORBIT": 11, "LWAVE": False, "NSW": 0, "ISYM": 0, "ICHARG": 11})
        self.kpoints_settings.update({"kpoints_density": kpoints_density})
        if mode == "Uniform":
            # Set smaller steps for DOS output
            self.incar_settings.update({"NEDOS": 601})
        if "NBANDS" not in user_incar_settings:
            raise KeyError("For NonSCF runs, NBANDS value from SC runs is "
                           "required!")
        else:
            self.incar_settings.update(user_incar_settings)

    def get_kpoints(self, structure):
        """
        Get a KPOINTS file for NonSCF calculation. In "Line" mode, kpoints are
        generated along high symmetry lines. In "Uniform" mode, kpoints are
        Gamma-centered mesh grid. Kpoints are written explicitly in both cases.

        Args:
            structure (Structure/IStructure): structure to get Kpoints
        """
        if self.mode == "Line":
            kpath = HighSymmKpath(structure)
            frac_k_points, k_points_labels = kpath.get_kpoints(
                line_density=self.kpoints_line_density,
                coords_are_cartesian=False)
            return Kpoints(comment="Non SCF run along symmetry lines",
                           style=Kpoints.supported_modes.Reciprocal,
                           num_kpts=len(frac_k_points),
                           kpts=frac_k_points, labels=k_points_labels,
                           kpts_weights=[1] * len(frac_k_points))
        else:
            num_kpoints = self.kpoints_settings["kpoints_density"] * \
                structure.lattice.reciprocal_lattice.volume
            kpoints = Kpoints.automatic_density(
                structure, num_kpoints * structure.num_sites)
            mesh = kpoints.kpts[0]
            ir_kpts = SpacegroupAnalyzer(structure, symprec=self.sym_prec) \
                .get_ir_reciprocal_mesh(mesh)
            kpts = []
            weights = []
            for k in ir_kpts:
                kpts.append(k[0])
                weights.append(int(k[1]))
            return Kpoints(comment="Non SCF run on uniform grid",
                           style=Kpoints.supported_modes.Reciprocal,
                           num_kpts=len(ir_kpts),
                           kpts=kpts, kpts_weights=weights)

    @staticmethod
    def get_incar_settings(vasp_run, outcar=None):
        """
        Helper method to get necessary user_incar_settings from previous run.

        Args:
            vasp_run (Vasprun): Vasprun that contains the final
                structure from previous run.
            outcar (Outcar): Outcar that contains the magnetization info
                from previous run.

        """
        # Turn off spin when magmom for every site is smaller than 0.02.
        if outcar and outcar.magnetization:
            site_magmom = np.array([i['tot'] for i in outcar.magnetization])
            ispin = 2 if np.any(site_magmom[np.abs(site_magmom) > 0.02]) else 1
        elif vasp_run.is_spin:
            ispin = 2
        else:
            ispin = 1
        nbands = int(np.ceil(vasp_run.as_dict()["input"]["parameters"]["NBANDS"]
                             * 1.2))
        incar_settings = {"ISPIN": ispin, "NBANDS": nbands}
        for grid in ["NGX", "NGY", "NGZ"]:
            if vasp_run.incar.get(grid):
                incar_settings.update({grid: vasp_run.incar.get(grid)})
        return incar_settings

    def get_incar(self, structure):
        incar = super(MPNonSCFVaspInputSet, self).get_incar(structure)
        incar.pop("MAGMOM", None)
        return incar

    def get_poscar(self, structure, get_primitive_standard=False):
        """
        Get a POSCAR file of the giving structure.

        Args:
            structure (Structure/IStructure): structure to get POSCAR
            get_primitive_standard (bool): if convert the input structure to a
            primitive standard structure
        """
        if get_primitive_standard:
            sym_finder = SpacegroupAnalyzer(structure, symprec=self.sym_prec)
            return Poscar(sym_finder.get_primitive_standard_structure(False))
        else:
            return Poscar(structure)

    @staticmethod
    def from_previous_vasp_run(previous_vasp_dir, output_dir='.',
                               mode="Uniform", user_incar_settings=None,
                               copy_chgcar=True, make_dir_if_not_present=True,
                               kpoints_density=1000, kpoints_line_density=20):
        """
        Generate a set of Vasp input files for NonSCF calculations from a
        directory of previous static Vasp run.

        Args:
            previous_vasp_dir (str): The directory contains the outputs(
                vasprun.xml and OUTCAR) of previous vasp run.
            output_dir (str): The directory to write the VASP input files
                for the NonSCF calculations. Default to write in the current
                directory.
            mode (str): Line: Generate k-points along symmetry lines for
                bandstructure. Uniform: Generate uniform k-points
                grids for DOS.
            user_incar_settings (dict): A dict specify customized settings
                for INCAR. Must contain a NBANDS value, suggest to use
                1.2*(NBANDS from static run).
            copy_chgcar (bool): Default to copy CHGCAR from SC run
            make_dir_if_not_present (bool): Set to True if you want the
                directory (and the whole path) to be created if it is not
                present.
            kpoints_density (int): kpoints density for the reciprocal cell
                of structure. Might need to increase the default value when
                calculating metallic materials.
            kpoints_line_density (int): kpoints density to use in line-mode.
                Might need to increase the default value when calculating
                metallic materials.
        """
        user_incar_settings = user_incar_settings or {}

        try:
            vasp_run = Vasprun(os.path.join(previous_vasp_dir, "vasprun.xml"),
                               parse_dos=False, parse_eigen=None)
            outcar = Outcar(os.path.join(previous_vasp_dir, "OUTCAR"))
            previous_incar = vasp_run.incar
        except:
            traceback.print_exc()
            raise RuntimeError("Can't get valid results from previous run: {}"
                               .format(previous_vasp_dir))

        #Get a Magmom-decorated structure
        structure = MPNonSCFVaspInputSet.get_structure(vasp_run, outcar,
                                                       initial_structure=True)
        nscf_incar_settings = MPNonSCFVaspInputSet.get_incar_settings(vasp_run,
                                                                      outcar)
        mpnscfvip = MPNonSCFVaspInputSet(nscf_incar_settings, mode,
                                         kpoints_density=kpoints_density,
                                         kpoints_line_density=kpoints_line_density)
        mpnscfvip.write_input(structure, output_dir, make_dir_if_not_present)
        if copy_chgcar:
            try:
                shutil.copyfile(os.path.join(previous_vasp_dir, "CHGCAR"),
                                os.path.join(output_dir, "CHGCAR"))
            except Exception as e:
                traceback.print_exc()
                raise RuntimeError("Can't copy CHGCAR from SC run" + '\n'
                                   + str(e))

        #Overwrite necessary INCAR parameters from previous runs
        previous_incar.update({"IBRION": -1, "ISMEAR": 0, "SIGMA": 0.001,
                               "LCHARG": False, "LORBIT": 11, "LWAVE": False,
                               "NSW": 0, "ISYM": 0, "ICHARG": 11})
        previous_incar.update(nscf_incar_settings)
        previous_incar.update(user_incar_settings)
        previous_incar.pop("MAGMOM", None)
        previous_incar.write_file(os.path.join(output_dir, "INCAR"))

        # Perform checking on INCAR parameters
        if any([previous_incar.get("NSW", 0) != 0,
                previous_incar["IBRION"] != -1,
                previous_incar["ICHARG"] != 11,
               any([sum(previous_incar["LDAUU"]) <= 0,
                    previous_incar["LMAXMIX"] < 4])
               if previous_incar.get("LDAU") else False]):
            raise ValueError("Incompatible INCAR parameters!")
