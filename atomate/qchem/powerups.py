# coding: utf-8


from atomate.qchem.firetasks.run_calc import RunQChemFake

__author__ = "Brandon Wood"
__email__ = "b.wood@berkeley.edu"


def use_fake_qchem(original_wf, ref_dirs, input_file="mol.qin"):
    """
        Replaces all RunQChem commands (i.e. RunQChemDirect, RunQChemCustodian) with RunQChemFake.
        This allows for testing without actually running QChem

        Args:
            original_wf (Workflow)
            ref_dirs (dict): key=firework name, value=path to the reference QChem calculation directory

        Returns:
            Workflow
    """
    for idx_fw, fw in enumerate(original_wf.fws):
        for job_type in ref_dirs.keys():
            if job_type == fw.name:
                for idx_t, t in enumerate(fw.tasks):
                    if "RunQChemCustodian" in str(t) or "RunQChemDirect" in str(t):
                        original_wf.fws[idx_fw].tasks[idx_t] = RunQChemFake(ref_dir=ref_dirs[job_type], input_file=input_file)

    return original_wf
