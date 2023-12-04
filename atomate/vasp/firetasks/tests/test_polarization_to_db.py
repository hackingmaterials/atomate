import os

from atomate.utils.testing import AtomateTest
from atomate.vasp.firetasks.parse_outputs import PolarizationToDb

__author__ = "Tess Smidt"
__email__ = "blondegeek@gmail.com"

module_dir = os.path.dirname(os.path.abspath(__file__))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
ref_dir = os.path.join(module_dir, "..", "..", "test_files")


DEBUG_MODE = (
    True  # If true, retains the database and output dirs at the end of the test
)
# If None, runs a "fake" VASP. Otherwise, runs VASP with this command...
VASP_CMD = None


class TestFerroelectricWorkflow(AtomateTest):
    def test_polarizationtodb(self):
        import gzip

        import bson

        reference_dir = os.path.abspath(os.path.join(ref_dir, "ferroelectric_wf"))

        with gzip.open(os.path.join(reference_dir, "tasks.bson.gz")) as f:
            coll_raw = f.read()

        coll = bson.decode_all(coll_raw)

        task_coll = self.get_task_collection()
        for c in coll:
            task_coll.insert_one(c)

        new_fw_spec = {
            "_fw_env": {"db_file": os.path.join(db_dir, "db.json")},
            "tags": ["wfid_1494203093.06934658"],
        }

        analysis = PolarizationToDb(db_file=">>db_file<<")
        analysis.run_task(new_fw_spec)

        # Check recovered change in polarization
        coll = self.get_task_collection("polarization_tasks")
        d = coll.find_one()
        self.assertAlmostEqual(d["polarization_change_norm"], 46.288752795325244, 5)
