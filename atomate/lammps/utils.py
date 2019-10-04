# coding: utf-8


from atomate.lammps.firetasks.run_calc import RunLammpsFake

__author__ = 'Kiran Mathew'
__email__ = "kmathew@lbl.gov"


def use_fake_lammps(original_wf, ref_dir):
    for idx_fw, fw in enumerate(original_wf.fws):
        for idx_t, t in enumerate(fw.tasks):
            if "RunLammps" in str(t):
                original_wf.fws[idx_fw].tasks[idx_t] = RunLammpsFake(ref_dir=ref_dir)
    return original_wf
