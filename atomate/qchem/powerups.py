from atomate.qchem.firetasks.run_calc import RunQChemFake

__author__ = "Brandon Wood", "Samuel Blau"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "11/20/19"


def use_fake_qchem(original_wf, ref_dirs, input_file="mol.qin"):
    """
    Replaces all RunQChem commands (i.e. RunQChemDirect, RunQChemCustodian) with RunQChemFake.
    This allows for testing without actually running QChem.
    Also deletes any RunCritic2 firetasks. Critic2 outputs will always be copied into the calc
    directory along with QChem outputs by RunQChemFake, so ProcessCritic2 will still run as expected.


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
                        original_wf.fws[idx_fw].tasks[idx_t] = RunQChemFake(
                            ref_dir=ref_dirs[job_type], input_file=input_file
                        )
                    if "RunCritic2" in str(t):
                        del original_wf.fws[idx_fw].tasks[idx_t]

    return original_wf
