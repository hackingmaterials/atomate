# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
This module defines a workflow for charged defects in non-metals:
- bulk calculation
- dielectric calculation (with option to omit)
- many different defects in several charge states (variety + charges used can be varied signficantly)

"""

from fireworks import Workflow, Firework

from pymatgen.io.vasp.sets import MPRelaxSet, MVLScanRelaxSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from atomate.vasp.powerups import add_modify_incar
from atomate.utils.utils import get_logger
from atomate.vasp.fireworks.core import OptimizeFW, StaticFW, HSEBSFW, DFPTFW
from atomate.vasp.firetasks.glue_tasks import CopyVaspOutputs

from atomate.vasp.firetasks.defects import DefectSetupFiretask, get_fw_from_defect

logger = get_logger(__name__)


def get_wf_chg_defects(structure, name="chg_defect_wf", user_incar_settings={},
                        vasp_cmd=">>vasp_cmd<<", db_file=">>db_file<<",
                        conventional=True, diel_flag=True, n_max=128,
                        vacancies=[], substitutions={}, interstitials={},
                        initial_charges={}, rerelax_flag=False, hybrid_run_for_gap_corr=True,
                        user_kpoints_settings={}):
    """
    Returns a charged defect workflow

    Firework 0 : (optional) re-relax the bulk structure that is input before running rest of workflow
    Firework 1 : (optional) run hybrid calculation of bulk structure to allow for
    Firework 2 : bulk supercell calculation
    Firework 3 : (optional) dielectric calculation
    Firework 4 - len(defectcalcs): Optimize the internal structure (fixed volume)
                                   of each charge+defect combination.

    (note if no re-relaxation required, then 1-len(defectcalcs) will all run at same time...)

    Args:
        structure (Structure): input structure to have defects run on
        name (str): some appropriate name for the workflow.
        user_incar_settings (dict):
            a dictionary of incar settings specified by user for both bulk and defect supercells
        vasp_cmd (str): Command to run vasp.
        db_file (str): path to file containing the database credentials.
        conventional (bool): flag to use conventional structure (rather than primitive) for defect supercells,
            defaults to True.
        diel_flag (bool): flag to also run dielectric calculations.
            (required for charge corrections to be run) defaults to True.
        n_max (int): maximum supercell size to consider for supercells

        vacancies (list):
            If list is totally empty, all vacancies are considered (default).
            If only specific vacancies are desired then add desired Element symbol to the list
                ex. ['Ga'] in GaAs structure will only produce Galium vacancies

            if NO vacancies are desired, then just add an empty list to the list
                ex. [ [] ]  yields no vacancies

        substitutions (dict):
            If dict is totally empty, all intrinsic antisites are considered (default).
            If only specific antisites/substituions are desired then add vacant site type as key, with list of
                sub site symbol as value
                    ex 1. {'Ga': ['As'] } in GaAs structure will only produce Arsenic_on_Gallium antisites
                    ex 2. {'Ga': ['Sb'] } in GaAs structure will only produce Antimonide_on_Gallium substitutions

            if NO antisites or substitutions are desired, then just add an empty dict
                ex. {'None':{}}  yields no antisites or subs


        interstitials (list):
            If list is totally empty, NO interstitial defects are considered (default).
            Option 1 for generation: If one wants to use Pymatgen to predict interstitial
                    then list of pairs of [symbol, generation method (str)] can be provided
                        ex. ['Ga', 'Voronoi'] in GaAs structure will produce Galium interstitials from the
                            Voronoi site finding algorithm
                        NOTE: only options for interstitial generation are "Voronoi" and "Nils"
            Option 2 for generation: If user wants to add their own interstitial sites for consideration
                    the list of pairs of [symbol, Interstitial object] can be provided, where the
                    Interstitial pymatgen.analysis.defects.core object is used to describe the defect of interest


        initial_charges (dict):
            says how to specify initial charges for each defect.
            An empty dict (DEFAULT) is to do a fairly restrictive charge generation method:
                for vacancies: use bond valence method to assign oxidation states and consider
                    negative of the vacant site's oxidation state as single charge to try
                antisites and subs: use bond valence method to assign oxidation states and consider
                    negative of the vacant site's oxidation state as single charge to try +
                    added to likely charge of substitutional site (closest to zero)
                interstitial: charge zero
            For non empty dict, charges are specified as:
                initial_charges = {'vacancies': {'Ga': [-3,2,1,0]},
                                   'substitutions': {'Ga': {'As': [0]} },
                                   'interstitials': {}}
                in the GaAs structure this makes vacancy charges in states -3,-2,-1,0; Ga_As antisites in the q=0 state,
                and all other defects will have charges generated in the restrictive automated format stated for DEFAULT

        rerelax_flag (bool):
            Flag to re-relax the input structure for minimizing forces 
            (does volume relaxation of small primitive cell) 
            Default is False (no re-relaxation occurs)
        hybrid_run_for_gap_corr (bool):
            Flag to run a single small hybrid bulk structure for band edge shifting correction.
            Default is True (Hybrid will be calculated)
        user_kpoints_settings (dict or Kpoints pmg object):
            a dictionary of kpoint settings specific by user OR an Actual Kpoint set to be used for the supercell calculations

    Returns:
        Workflow
    """
    fws, parents = [], []

    #force optimization and dielectric calculations with primitive structure for expediency
    prim_structure = SpacegroupAnalyzer(structure).find_primitive()

    if rerelax_flag:
        incar_settings = {"ISIF": 3, "ISPIN":2, "ISYM":2, "LAECHG":False, 'LVHAR': False}
        kpt_density = 100
        vis = MPRelaxSet( prim_structure,
                          user_incar_settings=incar_settings,
                          user_kpoints_settings={"reciprocal_density": kpt_density})

        rerelax_fw = OptimizeFW( prim_structure,
                                 name="structure optimization",
                                 vasp_input_set=vis,
                                 vasp_cmd=vasp_cmd, db_file=db_file,
                                 auto_npar=">>auto_npar<<",
                                 parents=None)
        fws.append(rerelax_fw)
        parents = [rerelax_fw]

        if hybrid_run_for_gap_corr: #only run HSEBSFW hybrid workflow here if re-relaxed since it requires a copy-over optimized structure
            hse_fw = HSEBSFW(structure=prim_structure, parents=parents, name="hse-BS", vasp_cmd=vasp_cmd, db_file=db_file)
            fws.append( hse_fw)

    elif hybrid_run_for_gap_corr: #if not re-relaxing structure but want hybrid then need to run a static primitive struct calc initial
        stat_gap_fw = StaticFW(structure=prim_structure, name="{} gap gga initialize".format(structure.composition.reduced_formula),
                                vasp_cmd=vasp_cmd, db_file=db_file)
        fws.append( stat_gap_fw)
        hse_fw = HSEBSFW(structure=prim_structure, parents=stat_gap_fw, name="hse-BS", vasp_cmd=vasp_cmd, db_file=db_file)
        fws.append( hse_fw)

    if diel_flag: #note dielectric DFPT run is only done with GGA
        user_incar_settings = {}
        copy_out = True if parents else False

        diel_fw = DFPTFW(structure=prim_structure, name='ionic dielectric', vasp_cmd=vasp_cmd, copy_vasp_outputs=copy_out,
                         db_file=db_file, parents=parents, user_incar_settings=user_incar_settings)
        fws.append( diel_fw)

    t = []
    if parents:
        t.append(CopyVaspOutputs(calc_loc= True ))

    t.append(DefectSetupFiretask(structure=prim_structure, cellmax=n_max, conventional=conventional,
                                 vasp_cmd=vasp_cmd, db_file=db_file, user_incar_settings=user_incar_settings,
                                 user_kpoints_settings=user_kpoints_settings,
                                 vacancies=vacancies, substitutions=substitutions,
                                 interstitials=interstitials, initial_charges=initial_charges))

    fw_name = "{} Defect Supercell Setup".format(structure.composition.reduced_formula)
    setup_fw = Firework(t,parents = parents, name=fw_name)
    fws.append(setup_fw)

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)
    final_wf = Workflow(fws, name=wfname)

    return final_wf


def run_defect_resubmissions( dpd, name="chg_defect_wf", consider_charges=True,
                              consider_supercells=False, user_incar_settings={},
                              vasp_cmd=">>vasp_cmd<<", db_file=">>db_file<<",
                              n_max=400, hybrid_run_for_gap_corr=False,
                              user_kpoints_settings={}):
    """
    This takes a DefectPhaseDiagram object and resubmits additional Defect calculations as needed.
    Also useful for resubmitting a hybrid calculation and band edge considerations.

    It can do this in two ways:
        (a) Extra Charge submissions: run any charges suggested for followup by DefectPhaseDiagram
        (b) Larger supercell submissions: run larger supercells for certain charges as suggested by DefectPhaseDiagram

    Args:
        dpd (DefectPhaseDiagram): phase diagram to consider for analysis of resubmission jobs
        name (str): some appropriate name for the workflow.
        consider_charges (bool): Whether to consider rerunning additional charge states
            Defaults to True (fireworks generated for additional charge states)
        consider_supercells (bool): Whether to consider submitting additional supercell sizes
            Note that when this is True, if there is no bulk_supercell that matches the supercell sizes
            being submitted, a bulk supercell calculation is also submitted.
            Defaults to False (will not generate fireworks for additional supercell sizes)
        user_incar_settings (dict):
            a dictionary of incar settings specified by user for both bulk and defect supercells
        vasp_cmd (str): Command to run vasp.
        db_file (str): path to file containing the database credentials.
        n_max (int): maximum supercell size to consider for consider_supercells supercell routine
        hybrid_run_for_gap_corr (bool):
            Flag to run a single small hybrid bulk structure for band edge shifting correction.
            Useful if you forgot to run the hybrid in your first
            Default is False (Hybrid will be calculated)
        user_kpoints_settings (dict or Kpoints pmg object):
            a dictionary of kpoint settings specific by user OR an Actual Kpoint set to be used for the
            supercell calculations

    Returns:
        Workflow
    :return:
    """
    if not consider_charges and not consider_supercells and not hybrid_run_for_gap_corr:
        raise ValueError("Ony makes since to run get_defect_resubmissions for resubmissions of "
                         "charges, supercells or hybrids.")

    fws = []

    if consider_charges:
        rec_dict = dpd.suggest_charges()
        for defname, charge_list in rec_dict.items():
            defect_template = dpd.finished_charges[ defname][0].defect.copy()
            if "scaling_matrix" in dpd.finished_charges[defname][0].parameters:
                supercell_size = dpd.finished_charges[defname][0].parameters['scaling_matrix']
            else:
                raise ValueError("No scaling_matrix available in {} defect entry! "
                                 "Cannot submit additional charged defects.".format( defname))

            for charge in charge_list:
                defect = defect_template.copy()
                defect.set_charge( charge)

                fw = get_fw_from_defect( defect, supercell_size,
                                         user_kpoints_settings = user_kpoints_settings,
                                         user_incar_settings = user_incar_settings,
                                         db_file=db_file, vasp_cmd=vasp_cmd)
                fws.append( fw)

    if consider_supercells:
        rec_dict = dpd.suggest_larger_supercells()
        for defname, charge_list in rec_dict.items():
            defect_template = dpd.finished_charges[ defname][0].defect.copy()

            #TODO -> consider next sized supercell size..
            supercell_size = ?

            for charge in charge_list:
                defect = defect_template.copy()
                defect.set_charge( charge)

                fw = get_fw_from_defect( defect, supercell_size,
                                         user_kpoints_settings = user_kpoints_settings,
                                         user_incar_settings = user_incar_settings,
                                         db_file=db_file, vasp_cmd=vasp_cmd)
                fws.append( fw)

    if not len(fws):
        print("No fireworks generated from defect resubmission scheme")
        return None
    else:
        bs_struct = dpd.entries[0].defect.bulk_structure.copy()
        wfname = "{}:{}".format(bs_struct.composition.reduced_formula, name)
        final_wf = Workflow(fws, name=wfname)

        return final_wf