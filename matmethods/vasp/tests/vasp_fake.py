import os
import shutil

from fireworks import FireTaskBase, explicit_serialize, Workflow
from matmethods.vasp.firetasks.run_calc import RunVaspCustodian
from pymatgen.io.vasp import Incar, Kpoints, Poscar, Potcar

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
reference_dir = os.path.join(module_dir, "reference_files")
fake_dirs={"structure optimization": os.path.join(reference_dir, "Si_structure_optimization"), "static": os.path.join(reference_dir, "Si_static"),
           "nscf uniform": os.path.join(reference_dir, "Si_nscf_uniform"), "nscf line": os.path.join(reference_dir, "Si_nscf_line")}

@explicit_serialize
class RunVaspFake(FireTaskBase):

    required_params = ["fake_dir"]

    def run_task(self, fw_spec):
        self._verify_inputs()
        self._clear_inputs()
        self._generate_outputs()

    def _verify_inputs(self):
        user_incar = Incar.from_file(os.path.join(os.getcwd(), "INCAR"))
        ref_incar = Incar.from_file(os.path.join(self["fake_dir"], "inputs", "INCAR"))

        # perform some BASIC tests

        # check INCAR
        params_to_check = ["ISPIN", "ENCUT", "ISMEAR", "SIGMA", "IBRION"]
        defaults = {"ISPIN": 1, "ISMEAR": 1, "SIGMA": 2}
        for p in params_to_check:
            if user_incar.get(p, defaults.get(p)) != ref_incar.get(p, defaults.get(p)):
                raise ValueError("INCAR value of {} is inconsistent!".format(p))

        # check KPOINTS
        user_kpoints = Kpoints.from_file(os.path.join(os.getcwd(), "KPOINTS"))
        ref_kpoints = Kpoints.from_file(os.path.join(self["fake_dir"], "inputs", "KPOINTS"))
        if user_kpoints._style != ref_kpoints.style or user_kpoints.num_kpts != ref_kpoints.num_kpts:
            raise ValueError("KPOINT files are inconsistent!".format(p))

        # check POSCAR
        user_poscar = Poscar.from_file(os.path.join(os.getcwd(), "POSCAR"))
        ref_poscar = Poscar.from_file(os.path.join(self["fake_dir"], "inputs", "POSCAR"))
        if user_poscar.natoms != ref_poscar.natoms or user_poscar.site_symbols != ref_poscar.site_symbols:
            raise ValueError("POSCAR files are inconsistent!".format(p))

        # check POTCAR
        user_potcar = Potcar.from_file(os.path.join(os.getcwd(), "POTCAR"))
        ref_potcar = Potcar.from_file(os.path.join(self["fake_dir"], "inputs", "POTCAR"))
        if user_potcar.symbols != ref_potcar.symbols:
            raise ValueError("POTCAR files are inconsistent!".format(p))

        print("RunVaspFake: verified inputs successfully")

    def _clear_inputs(self):
        for x in ["INCAR", "KPOINTS", "POSCAR", "POTCAR", "CHGCAR", "OUTCAR", "vasprun.xml"]:
            p = os.path.join(os.getcwd(), x)
            if os.path.exists(p):
                os.remove(p)

    def _generate_outputs(self):
        output_dir = os.path.join(self["fake_dir"], "outputs")
        for file_name in os.listdir(output_dir):
            full_file_name = os.path.join(output_dir, file_name)
            if (os.path.isfile(full_file_name)):
                shutil.copy(full_file_name, os.getcwd())

        print("RunVaspFake: ran fake VASP, generated outputs")


def make_fake_workflow(original_workflow):
    wf_dict = original_workflow.to_dict()
    for idx_fw, fw in enumerate(original_workflow.fws):
        for job_type in fake_dirs.keys():
            if job_type in fw.name:
                for idx_t, t in enumerate(fw.tasks):
                    if "RunVasp" in str(t):
                        wf_dict["fws"][idx_fw]["spec"]["_tasks"][idx_t] = RunVaspFake(fake_dir=fake_dirs[job_type]).to_dict()

    return Workflow.from_dict(wf_dict)

# TODO: move this to "powerups"
def make_custodian_workflow(original_workflow):
    wf_dict = original_workflow.to_dict()
    for idx_fw, fw in enumerate(original_workflow.fws):
        for job_type in fake_dirs.keys():
            if job_type in fw.name:
                for idx_t, t in enumerate(fw.tasks):
                    if "RunVasp" in str(t):
                        vasp_cmd = wf_dict["fws"][idx_fw]["spec"]["_tasks"][idx_t]["vasp_cmd"]
                        wf_dict["fws"][idx_fw]["spec"]["_tasks"][idx_t] = RunVaspCustodian(vasp_cmd=vasp_cmd).to_dict()

    return Workflow.from_dict(wf_dict)