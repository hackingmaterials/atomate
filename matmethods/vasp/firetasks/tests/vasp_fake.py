import os
import shutil

from fireworks import FireTaskBase, explicit_serialize

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


@explicit_serialize
class VaspFakeTask(FireTaskBase):

    required_params = ["run_dir"]

    def run_task(self, fw_spec):
        self._verify_inputs()
        self._generate_outputs()


    def _verify_inputs(self):
        pass
        print("RAN VERIFY INPUTS")

    def _generate_outputs(self):
        output_dir = os.path.join(self["run_dir"], "outputs")
        for file_name in os.listdir(output_dir):
            full_file_name = os.path.join(output_dir, file_name)
            if (os.path.isfile(full_file_name)):
                shutil.copy(full_file_name, os.getcwd())

        print("RAN GENERATE OUTPUTS")