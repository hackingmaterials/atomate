import os
import shutil

from fireworks import FireTaskBase, explicit_serialize
from pymatgen.io.vasp import Incar, Kpoints, Poscar, Potcar

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


@explicit_serialize
class RunVaspFake(FireTaskBase):

    required_params = ["run_dir"]

    def run_task(self, fw_spec):
        self._verify_inputs()
        self._generate_outputs()


    def _verify_inputs(self):
        user_incar = Incar.from_file(os.path.join(os.getcwd(), "INCAR"))
        ref_incar = Incar.from_file(os.path.join(self["run_dir"], "inputs", "INCAR"))

        # perform some BASIC tests

        # check INCAR
        params_to_check = ["ISPIN", "ENCUT", "ISMEAR", "SIGMA", "IBRION"]
        defaults = {"ISPIN": 1, "ISMEAR": 1, "SIGMA": 2}
        for p in params_to_check:
            if user_incar.get(p, defaults.get(p)) != ref_incar.get(p, defaults.get(p)):
                raise ValueError("INCAR value of {} is inconsistent!".format(p))

        # check KPOINTS
        user_kpoints = Kpoints.from_file(os.path.join(os.getcwd(), "KPOINTS"))
        ref_kpoints = Kpoints.from_file(os.path.join(self["run_dir"], "inputs", "KPOINTS"))
        if user_kpoints._style != ref_kpoints.style or user_kpoints.num_kpts != ref_kpoints.num_kpts:
            raise ValueError("KPOINT files are inconsistent!".format(p))

        # check POSCAR
        user_poscar = Poscar.from_file(os.path.join(os.getcwd(), "POSCAR"))
        ref_poscar = Poscar.from_file(os.path.join(self["run_dir"], "inputs", "POSCAR"))
        if user_poscar.natoms != ref_poscar.natoms or user_poscar.site_symbols != ref_poscar.site_symbols:
            raise ValueError("POSCAR files are inconsistent!".format(p))

        # check POTCAR
        user_potcar = Potcar.from_file(os.path.join(os.getcwd(), "POTCAR"))
        ref_potcar = Potcar.from_file(os.path.join(self["run_dir"], "inputs", "POTCAR"))
        if user_potcar.symbols != ref_potcar.symbols:
            raise ValueError("POTCAR files are inconsistent!".format(p))

        print("RunVaspFake: verified inputs successfully")


    def _generate_outputs(self):
        output_dir = os.path.join(self["run_dir"], "outputs")
        for file_name in os.listdir(output_dir):
            full_file_name = os.path.join(output_dir, file_name)
            if (os.path.isfile(full_file_name)):
                shutil.copy(full_file_name, os.getcwd())

        print("RunVaspFake: ran fake VASP, generated outputs")