# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import unittest

from fireworks import Firework, ScriptTask, Workflow

from atomate.vasp.powerups import add_priority, use_custodian, add_trackers, \
    add_modify_incar, add_small_gap_multiply, use_scratch_dir, remove_custodian, \
    add_tags, add_wf_metadata
from atomate.vasp.workflows.base.core import get_wf

from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.util.testing import PymatgenTest


__author__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'


class TestVaspPowerups(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        struct_si = PymatgenTest.get_structure("Si")
        vis = MPRelaxSet(struct_si, force_gamma=True)
        cls.bs_wf = get_wf(struct_si,
                           "bandstructure.yaml",
                           vis=vis, common_params={"vasp_cmd": "test_VASP"})
        cls.bsboltz_wf = get_wf(struct_si,
                           "bandstructure_boltztrap.yaml",
                           vis=vis)

    def _copy_wf(self, wf):
        return Workflow.from_dict(wf.to_dict())

    def test_add_priority(self):
        fw1 = Firework([ScriptTask(script=None)], fw_id=-1)
        fw2 = Firework([ScriptTask(script=None)], parents=[fw1], fw_id=-2)
        fw3 = Firework([ScriptTask(script=None)], parents=[fw1], fw_id=-3)

        wf = Workflow([fw1, fw2, fw3])

        wf = add_priority(wf, 4, 8)
        self.assertEqual(wf.id_fw[-1].spec["_priority"], 4)
        self.assertEqual(wf.id_fw[-2].spec["_priority"], 8)
        self.assertEqual(wf.id_fw[-3].spec["_priority"], 8)

    def test_custodian_powerups(self):
        my_wf = self._copy_wf(self.bs_wf)
        my_wf = remove_custodian(my_wf)

        for fw in my_wf.fws:
            task_idx = 1 if "structure optimization" in fw.name else 2
            self.assertTrue(
                "RunVaspDirect" in fw.tasks[task_idx]._fw_name)
            self.assertEqual(
                fw.tasks[task_idx]["vasp_cmd"], "test_VASP")

        my_wf_double_relax = remove_custodian(self._copy_wf(self.bs_wf))
        my_wf_double_relax = use_custodian(my_wf_double_relax,
                                           fw_name_constraint="structure optimization",
                                           custodian_params={"job_type": "double_relaxation_run"})

        for fw in my_wf_double_relax.fws:
            if "structure optimization" in fw.name:
                self.assertTrue("RunVaspCustodian" in fw.tasks[1]._fw_name)
                self.assertEqual(fw.tasks[1]["job_type"],
                                 "double_relaxation_run")
            else:
                self.assertTrue("RunVaspDirect" in fw.tasks[2]._fw_name)
                self.assertFalse("job_type" in fw.tasks[2])

    def test_modify_incar(self):
        my_wf = add_modify_incar(self._copy_wf(self.bs_wf), {"incar_update": {"NCORE": 1}},
                                 fw_name_constraint="structure optimization")

        for fw in my_wf.fws:
            if "structure optimization" in fw.name:
                self.assertTrue("ModifyIncar" in fw.tasks[1]._fw_name)
                self.assertEqual(fw.tasks[1]["incar_update"], {"NCORE": 1})
            else:
                for t in fw.tasks:
                    self.assertFalse("ModifyIncar" in t["_fw_name"])

    def test_add_trackers(self):
        my_wf = add_trackers(self._copy_wf(self.bs_wf))

        for fw in my_wf.fws:
            self.assertEqual(len(fw.spec["_trackers"]), 2)

    def test_add_small_gap_multiply(self):
        my_wf = self._copy_wf(self.bs_wf)
        my_wf = add_small_gap_multiply(my_wf, 0.5, 1.5, "static")
        found = False

        for fw in my_wf.fws:
            if "static" in fw.name:
                for t in fw.tasks:
                    if 'WriteVasp' in str(t):
                        self.assertEqual(t["small_gap_multiply"], [0.5, 1.5])
                        found = True

        self.assertEqual(found, True)

    def test_use_scratch_dir(self):
        my_wf = self._copy_wf(self.bs_wf)
        my_wf = use_custodian(my_wf)
        my_wf = use_scratch_dir(my_wf, ">>scratch_dir<<")
        found = 0

        for fw in my_wf.fws:
            for t in fw.tasks:
                if 'RunVaspCustodian' in str(t):
                    self.assertEqual(t["scratch_dir"], ">>scratch_dir<<")
                    found += 1

        self.assertEqual(found, 4)

    def test_add_tags(self):
        my_wf = self._copy_wf(self.bs_wf)
        my_wf.metadata = {"tags": ["a"]}
        my_wf = add_tags(my_wf, ["b", "c"])

        found = 0

        self.assertEqual(my_wf.metadata["tags"], ["a", "b", "c"])
        for fw in my_wf.fws:
            self.assertEqual(fw.spec["tags"], ["b", "c"])
            for t in fw.tasks:
                if 'VaspToDb' in str(t):
                    self.assertEqual(t["additional_fields"]["tags"], ["b", "c"])
                    found += 1
        self.assertEqual(found, 4)

        my_wf = self._copy_wf(self.bsboltz_wf)
        my_wf = add_tags(my_wf, ["foo", "bar"])

        v_found = 0
        b_found = 0

        self.assertEqual(my_wf.metadata["tags"], ["foo", "bar"])
        for fw in my_wf.fws:
            self.assertEqual(fw.spec["tags"], ["foo", "bar"])
            for t in fw.tasks:
                if 'BoltztrapToDb' in str(t):
                    self.assertEqual(t["additional_fields"]["tags"], ["foo", "bar"])
                    b_found += 1
                if 'VaspToDb' in str(t):
                    self.assertEqual(t["additional_fields"]["tags"], ["foo", "bar"])
                    v_found += 1
        self.assertEqual(b_found, 1)
        self.assertEqual(v_found, 4)

    def test_add_wf_metadata(self):
        my_wf = self._copy_wf(self.bs_wf)
        my_wf = add_wf_metadata(my_wf, PymatgenTest.get_structure("Si"))
        self.assertEqual(my_wf.metadata["nelements"], 1)
        self.assertEqual(my_wf.metadata["formula"], "Si2")


if __name__ == "__main__":
    unittest.main()
