import bson
import gzip

from pathlib import Path

from atomate.utils.testing import AtomateTest
from atomate.vasp.firetasks.parse_outputs import PolarizationToDb

__author__ = "Tess Smidt"
__email__ = "blondegeek@gmail.com"

module_dir = Path(__file__).resolve().parent
db_dir = module_dir / "../../../common/test_files"
ref_dir = module_dir / "../../test_files"


class TestPolarizationFiretasks(AtomateTest):
    def test_polarizationtodb(self):
        wf_dir = ref_dir / "ferroelectric_wf"

        with gzip.open(wf_dir / "tasks.bson.gz") as f:
            coll_raw = f.read()

        coll = bson.decode_all(coll_raw)

        db = self.get_task_collection()
        for c in coll:
            db.insert(c)

        new_fw_spec = {
            "_fw_env": {"db_file": db_dir / "db.json"},
            "tags": ["wfid_1494203093.06934658"],
        }

        analysis = PolarizationToDb(db_file=">>db_file<<")
        analysis.run_task(new_fw_spec)

        # Check recovered change in polarization
        coll = self.get_task_collection("polarization_tasks")
        d = coll.find_one()
        self.assertAlmostEqual(d["polarization_change_norm"], 46.28875279532, 5)
