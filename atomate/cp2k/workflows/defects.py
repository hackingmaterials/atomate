"""
This module defines a workflow for charged defects in non-metals. It is based on the original workflow
developed by Danny Broberg for VASP.
"""

from fireworks import Workflow, Firework
from copy import deepcopy
from pymatgen.io.vasp.sets import MPRelaxSet, MVLScanRelaxSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from atomate.utils.utils import get_logger
from atomate.cp2k.fireworks.core import StaticFW, RelaxFW, StaticHybridFW, RelaxHybridFW
from atomate.cp2k.utils import get_defect_structures, optimize_structure_sc_scale

#from atomate.vasp.firetasks.defects import DefectSetupFiretask
#from atomate.vasp.fireworks.defects import DefectAnalysisFW

logger = get_logger(__name__)


# TODO: does conventional unit cell vs. primitive matter?
# TODO: run_analysis FW for all defects in WF or should it be a part of DefectDrone
def get_wf_chg_defects(structure,
                       mpid=None,
                       name="chg_defect_wf",
                       cp2k_cmd=">>cp2k_cmd<<",
                       db_file=">>db_file<<",
                       cp2k_gga_input_set=None,
                       user_gga_settings={},
                       cp2k_hybrid_input_set=None,
                       user_hybrid_settings={},
                       cp2ktodb_kwargs={},
                       diel_flag=False,
                       defect_dict={},
                       rerelax_flag=False,
                       scale_to_num_sites=500):
    """
    Returns a charged defect workflow.

    Firework 0 : (optional) re-relax the bulk structure that is input before running rest of workflow
    Firework 1 : (optional) run hybrid calculation of bulk structure to allow for
    Firework 2 : bulk supercell calculation
    Firework 3 : (optional) dielectric calculation
    Firework 4 - len(defectcalcs): Optimize the internal structure (fixed volume)
                                   of each charge+defect combination.

    (note if no re-relaxation required, then 1-len(defectcalcs) will all run at same time...)

    Args:
        structure (Structure): input structure to have defects run on
        mpid (str): Materials Project id number to be used (for storage in metadata).
        name (str): some appropriate name for the workflow.
        cp2k_gga_input_set (Cp2kInputSet): User defined Cp2kInputSet to use for the GGA calculations
            in this workflow in place of the default.
            Default: None
        user_gga_settings (dict): Dictionary defining user settings for initializing the Cp2kInputSet
            for GGA calculations. See Cp2kInputSet for list of arguments given to Cp2kInputSet.
            This can be  more convenient to use this instead of cp2k_gga_input_set.
            Default: {}
        cp2k_hybrid_input_set (Cp2kInputSet): User defined Cp2kInputSet to use for the hybrid calculations
            in this workflow in place of the default.
            Default: None
        user_hybrid_settings (dict): Dictionary defining user settings for initializing the Cp2kInputSet
            for hybrid calculations. See Cp2kInputSet for list of arguments given to Cp2kInputSet.
            This can be  more convenient to use this instead of cp2k_gga_input_set.
            Default: {}
        cp2k_cmd (str): Command to run cp2k.
            Default: >>cp2k_cmd<< (use env check on fworker)
        db_file (str): path to file containing the database credentials.
            Default: >>db_file<< (use env check on fworker)
        conventional (bool): flag to use conventional structure (rather than primitive) for defect supercells,
            defaults to True.
        diel_flag (bool): flag to also run dielectric calculations.
            (required for charge corrections to be run) defaults to True.
        n_max (int): maximum supercell size to consider for supercells

        job_type (str): type of defect calculation that user desires to run
            default is 'normal' which runs a GGA defect calculation
            additional options are:
                'double_relaxation_run' which runs a double relaxation GGA run
                'metagga_opt_run' which runs a double relaxation with SCAN (currently
                    no way to turn off double relaxation approach with SCAN)
                'hse' which runs a relaxation step with GGA followed by a relaxation with HSE
                    NOTE: for these latter two we highly recommend that rerelax is set to True
                        so the bulk_structure's lattice is optimized before running defects

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
    Returns:
        Workflow
    """
    fws, parents = [], []

    # Make supercell. Should always be done for CP2K because it is gamma point only
    structure.make_supercell(optimize_structure_sc_scale(structure, scale_to_num_sites))

    # Use the DefectDrone unless specified
    cp2ktodb_kwargs['drone'] = cp2ktodb_kwargs.get('drone', 'DefectDrone')

    # Re-relax the initial structure
    if rerelax_flag:
        fws.append(
            RelaxFW(
                structure=structure,
                name='Re-Relax-FW',
                cp2k_input_set=cp2k_gga_input_set,
                cp2k_cmd=cp2k_cmd,
                prev_calc_loc=False,
                db_file=db_file,
                cp2ktodb_kwargs=cp2ktodb_kwargs,
                parents=None
            )
        )
        parents.append(fws[-1])

    # Run dielectric calculation before defects (for finite size corrections)
    if diel_flag:
        print("DIELECTRIC FW FOR CP2K IS NOT IMPLEMENTED")

    # Collect all defects based on the defect provided and run GGA-->Hybrid for each
    defects = get_defect_structures(structure, defect_dict=defect_dict)
    for i, defect in enumerate(defects):
        bulk_name = 'Re-Relax-FW' if rerelax_flag else None  # So prev_calc_loc can find re-relax fw
        cp2ktodb_kwargs['fw_spec_field'].update({'defect': defect.as_dict()})  # Keep record of defect object
        gga_name = "Defect-GGA-FW-{}".format(i)  # How to track the GGA FW
        hybrid_name = "Defect-Hybrid-FW-{}".format(i)  # How to track the hybrid FW
        fws.append(
            RelaxFW(
                structure=defect.bulk_structure,
                name=gga_name,
                cp2k_input_set=cp2k_gga_input_set,
                cp2k_input_set_params=user_gga_settings,
                cp2k_cmd=cp2k_cmd,
                prev_calc_loc=bulk_name,
                db_file=db_file,
                cp2ktodb_kwargs=cp2ktodb_kwargs,
                parents=parents,
                files_to_copy=None
            )
        )
        fws.append(
            StaticHybridFW(
                structure=defect.bulk_structure,
                name=hybrid_name,
                cp2k_input_set=cp2k_hybrid_input_set,
                cp2k_input_set_params=user_hybrid_settings,
                cp2k_cmd=cp2k_cmd,
                prev_calc_loc=gga_name,
                db_file=db_file,
                cp2ktodb_kwargs=cp2ktodb_kwargs,
                parents=fws[-1],
                files_to_copy="{}-RESTART.wfn".format(gga_name)
            )
        )

    wfname = "{}:{}".format(structure.composition.reduced_formula, name)
    final_wf = Workflow(fws, name=wfname)

    return final_wf


