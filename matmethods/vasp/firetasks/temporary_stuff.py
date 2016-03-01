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
    def write_input_from_prevrun():
        # TODO: implement me!
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

        """
        pass