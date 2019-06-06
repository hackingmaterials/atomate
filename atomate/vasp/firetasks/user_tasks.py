import os
import random

from pymatgen.io.lammps.data import LammpsData
from pymatgen.io.vasp.sets import MITMDSet

from fireworks.core.firework import FiretaskBase, FWAction
from fireworks import explicit_serialize

from atomate.utils.utils import get_logger
from atomate.vasp.fireworks.core import MDFW

logger = get_logger(__name__)

@explicit_serialize
class LammpsToVaspMD(FiretaskBase):
    _fw_name = "LammpsToVasp"
    required_params = ["atom_style", "start_temp", "end_temp", "nsteps"]
    optional_params = ['time_step', 'vasp_input_set', 'user_kpoints_settings', 'vasp_cmd',
                       'copy_vasp_outputs', 'db_file', 'name', 'parents']

    def run_task(self, fw_spec):
        atom_style  = self.get('atom_style')
        start_temp  = self.get('start_temp')
        end_temp    = self.get('end_temp')
        nsteps      = self.get('nsteps')

        time_step = self.get('time_step') or 1
        vasp_cmd = self.get('vasp_cmd') or ">>vasp_cmd<<"
        copy_vasp_outputs = self.get('copy_vasp_outputs') or False
        db_file = self.get('db_file') or None
        name = self.get('name') or "VaspMDFW"
        parents = self.get('parents') or None
        transmute = self.get('transmute') or None

        logger.info("PARSING \"lammps.final\" to VASP.")
        data = LammpsData.from_file(os.path.join(os.getcwd(), self.get('final_data')),
                                    atom_style=atom_style, sort_id=True)
        structure = data.structure

        if transmute:
            sites = structure.sites
            indices = []
            for i, s in enumerate(sites):
                if s.specie.symbol == transmute[0]:
                    indices.append(i)
            structure.replace(random.choice(indices), species=transmute[1])

        vasp_input_set = fw_spec.get('vasp_input_set') or MITMDSet(structure, start_temp, end_temp, nsteps, time_step,
                                                                   force_gamma=True)
        user_kpoints_settings = fw_spec.get('user_kpoints_settings') or None

        if user_kpoints_settings:
            v = vasp_input_set.as_dict()
            v.update({"user_kpoints_settings": user_kpoints_settings})
            vasp_input_set = vasp_input_set.from_dict(v)

        fw = MDFW(structure, start_temp, end_temp, nsteps, vasp_input_set=vasp_input_set, vasp_cmd=vasp_cmd,
                  copy_vasp_outputs=copy_vasp_outputs, db_file=db_file, name=name)

        return FWAction(additions=fw)