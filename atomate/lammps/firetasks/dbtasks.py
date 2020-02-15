import json
import os
import re
import gridfs
import zlib

import numpy as np
from bson import ObjectId
from atomate.utils.utils import env_chk, get_logger
from atomate.lammps.database import LammpsCalcDb
from fireworks import explicit_serialize, FiretaskBase
from monty.json import MontyEncoder
from pymatgen.core.trajectory import Trajectory

__author__ = 'Eric Sivonxay'

logger = get_logger(__name__)

@explicit_serialize
class LammpsMDToDB(FiretaskBase):
    """
    Obtain all production runs and insert them into the db. This is done by
    searching for a unique tag
    """
    required_params = ["db_file"]
    optional_params = ["input_filename"]

    def run_task(self, fw_spec):
        db_file = env_chk(self["db_file"], fw_spec)

        input_filename = self.get("input_filename", 'lammps.in')

        with open(input_filename, 'r') as f:
            input_string = f.read()

        timestep = float(re.search(r"timestep\s+(.*)\n", input_string).group(1).split()[0])
        temp = float(re.search(r"variable\s+ mytemp(.*)\n", input_string).group(1).split()[1])

        groups = re.findall(r"group\s+(.*)\n", input_string)
        atomic_map = {}
        for group_str in groups:
            atom, _str, atom_id = group_str.split()
            atomic_map[int(atom_id)] = atom

        # Load trajectories
        path = os.getcwd()
        dump_filename = re.search(r"variable\s+ NAME(.*)\n", input_string).group(1).split()[1]
        dumpfile = os.path.join(path, f'{dump_filename}.lammpstrj')
        trajectory = LammpsTrajectory.from_lammps_dump(dumpfile, atomic_map, time_step=timestep)

        mmdb = LammpsCalcDb.from_db_file(db_file, admin=True)
        traj_dict = json.dumps(trajectory.as_dict(), cls=MontyEncoder)
        gfs_id, compression_type = insert_gridfs(traj_dict, mmdb.db, "trajectories_fs")

        # insert trajectory into db

        traj_doc = {
            'formula_pretty': trajectory[0].composition.reduced_formula,
            'formula': trajectory[0].composition.formula.replace(' ', ''),
            'temperature': temp,
            'compression': compression_type,
            'fs_id': gfs_id,
            'structure': trajectory[0].as_dict(),
            'dimension': list(np.shape(trajectory.frac_coords)),
            'time_step': timestep
        }

        mmdb.db.trajectories.insert_one(traj_doc)


def insert_gridfs(d, db, collection="fs", compress=True, oid=None, task_id=None):
    """
    Insert the given document into GridFS.
    Args:
        d (dict): the document
        collection (string): the GridFS collection name
        compress (bool): Whether to compress the data or not
        oid (ObjectId()): the _id of the file; if specified, it must not already exist in GridFS
        task_id(int or str): the task_id to store into the gridfs metadata
    Returns:
        file id, the type of compression used.
    """
    oid = oid or ObjectId()
    compression_type = None

    if compress:
        d = zlib.compress(d.encode(), compress)
        compression_type = "zlib"

    fs = gridfs.GridFS(db, collection)
    if task_id:
        # Putting task id in the metadata subdocument as per mongo specs:
        # https://github.com/mongodb/specifications/blob/master/source/gridfs/gridfs-spec.rst#terms
        fs_id = fs.put(d, _id=oid, metadata={"task_id": task_id, "compression": compression_type})
    else:
        fs_id = fs.put(d, _id=oid, metadata={"compression": compression_type})

    return fs_id, compression_type


class LammpsTrajectory(Trajectory):
    #     def __init__(avg_msd, )

    @classmethod
    def from_lammps_dump(cls, file, species_map, time_step=1):
        with open(file, 'r') as f:
            lines = f.readlines()

        #### Assume number of atoms doesn't change
        section_len = None
        timestep_i = None
        natoms_i = None
        box_i = None
        atoms_i = None

        for i, line in enumerate(lines):
            # Find length of each section
            if line.startswith('ITEM: TIMESTEP'):
                if i > 1:
                    section_len = i
                    break
                else:
                    timestep_i = i + 1
            elif line.startswith('ITEM: NUMBER'):
                natoms_i = i + 1
            elif line.startswith('ITEM: BOX'):
                box_i = i + 1
            elif line.startswith('ITEM: ATOMS'):
                atoms_i = i + 1

        time_steps = []
        natoms = []
        box_dims = []
        coordinates = []
        species = []

        for i in range(0, len(lines), section_len):
            time_steps.append(re.split('\s+', lines[i + timestep_i].strip())[0])
            natom = int(re.split('\s+', lines[i + natoms_i].strip())[0])
            natoms.append(natom)

            # Parse box size
            box = []
            for l in lines[i + box_i: i + box_i + 3]:
                box.append([float(i) for i in re.split('\s+', l.strip())])
            lattice = [[box[0][1] - box[0][0], 0, 0], [0, box[1][1] - box[1][0], 0], [0, 0, box[2][1] - box[2][0]]]
            box_dims.append(lattice)

            # Parse coordinates
            _coordinates = []
            for l in lines[i + atoms_i: i + atoms_i + natom]:
                _coordinates.append([float(i) for i in re.split('\s+', l.strip())])
            spec_coords = np.array(sorted(_coordinates, key=lambda x: x[0]))
            coordinates.append(np.divide(spec_coords[:, 2:], np.diag(lattice)))

            if i == 0:
                species = list(map(int, spec_coords[:, 1]))
        time_steps = list(map(int, time_steps))
        species = [species_map[i] for i in species]
        return cls(lattice=box_dims, species=species,
                   frac_coords=coordinates, constant_lattice=False,
                   time_step=(time_steps[1] - time_steps[0]) * time_step)

