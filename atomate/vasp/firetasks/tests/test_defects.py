
import numpy as np

import unittest

from atomate.vasp.firetasks.defects import optimize_structure_sc_scale, DefectSetupFiretask

from pymatgen.core import Structure, PeriodicSite
from pymatgen.analysis.defects.core import Interstitial
from pymatgen.util.testing import PymatgenTest


class TestOptimize(PymatgenTest):

    def test_optimize_structure_sc_scale(self):
        #basic cubic structure test
        s = PymatgenTest.get_structure("CsCl")
        lattchange = optimize_structure_sc_scale(s, 300)
        self.assertSequenceEqual([5, 5, 5], lattchange)
        lattchange = optimize_structure_sc_scale(s, 100)
        self.assertSequenceEqual([3, 3, 3], lattchange)

        #non cubic test
        s = PymatgenTest.get_structure("Si")
        lattchange = optimize_structure_sc_scale(s, 200)
        self.assertSequenceEqual([4, 4, 4], lattchange)
        lattchange = optimize_structure_sc_scale(s, 100)
        self.assertSequenceEqual([3, 3, 3], lattchange)

        #test asking for smaller supercell size than original structure
        s.make_supercell(2)
        lattchange = optimize_structure_sc_scale(s, len(s) - 2)
        self.assertSequenceEqual([1, 1, 1], lattchange)

        #test additional scaling for anisotropic lattice constants
        s = PymatgenTest.get_structure("Graphite")
        lattchange = optimize_structure_sc_scale(s, 300)
        self.assertSequenceEqual([6, 6, 2], lattchange)
        lattchange = optimize_structure_sc_scale(s, 100)
        self.assertSequenceEqual([3, 3, 2], lattchange)

        #test additional scaling for anisotropic lattice angles
        s = PymatgenTest.get_structure("Li2O")
        lattchange = optimize_structure_sc_scale(s, 300)
        self.assertSequenceEqual([4, 4, 4], lattchange)
        lattchange = optimize_structure_sc_scale(s, 100)
        self.assertSequenceEqual([3, 3, 3], lattchange)


class TestDefectSetupFiretask(PymatgenTest):

    def _verify_setuptask(self, setuptask, supercell=[1.,1.,1.], user_incar_settings={},
                          charge=0, is_defect=True):
        """
        verify that charges, incar settings, supercell scaling,
        and charge were set correctly for each setup task
        """
        self.assertArrayEqual( setuptask['transformation_params'][0]['scaling_matrix'], supercell)

        #default incar settings used for bulk and defect calcs
        incar = {"EDIFF":.0001, "EDIFFG": 0.001, "ISMEAR":0, "SIGMA":0.05, "ISIF": 2,
                         "ISPIN":2,  "ISYM":2, "LVHAR":True, "LVTOT":True, "LWAVE": True}

        if is_defect:
            incar.update( {"NSW": 100, "LAECHG":False})
            self.assertEqual( setuptask['vasp_input_set']['structure']['charge'], charge)
            self.assertTrue( setuptask['vasp_input_set']['use_structure_charge'])
        else:
            incar.update( {"NSW": 0})

        incar.update(user_incar_settings)
        for check_key, check_val in incar.items():
            self.assertEqual( setuptask['vasp_input_set']['user_incar_settings'][check_key], check_val)


    def _verify_defect_list(self, defect_setuptask_list, defects_expected_list):
        """
        verify that defect list that was setup is equal to expected defect list.

        :param defect_setuptask_list:
        :param defects_expected_list: organized as a list of [DefectClassString, Specie_to_add, Specie_to_remove, charge]
            examples: a neutral Vacancy = ['Vacancy', None, 'Cs', 0]
                      a neutral Substitution of Cs_on_Cl is ['Substitution', 'Cs', 'Cl', 0]
                      a neutral Interstitial with Cs is ['Interstitial', None, 'Cl', 0]
        :return:
        """
        #copy expected defect list and remove items from this list
        # as they are found in the defect_setuptask_list. Fails if list is not empty
        cop_def_exp = defects_expected_list[:]

        #each defect setup task has keys = ['@class', 'structure', 'charge', 'defect_site']
        for d in defect_setuptask_list:
            charge = d['charge']
            sub_site = d['defect_site']['species'][0]['element']
            if d['@class'] == 'Vacancy':
                this_def = ['Vacancy', None, sub_site, charge]
            elif d['@class'] == 'Interstitial':
                this_def = ['Interstitial', sub_site, None, charge]
            elif d['@class'] == 'Substitution':
                tmp_struct = Structure.from_dict( d['structure'])
                tmp_site = PeriodicSite.from_dict( d['defect_site'])
                poss_deflist = sorted(
                    tmp_struct.get_sites_in_sphere(tmp_site.coords, 2, include_index=True), key=lambda x: x[1])
                defindex = poss_deflist[0][2]
                vac_site = tmp_struct[defindex].specie.symbol
                this_def = ['Substitution', sub_site, vac_site, charge]

            cop_def_exp.remove( this_def)

        self.assertEqual( len(cop_def_exp), 0)


    def test_defect_setup(self):
        struct = PymatgenTest.get_structure("CsCl")

        #test standard defect setup with defaults
        ft = DefectSetupFiretask( structure=struct, cellmax=100)

        def_set = ft.run_task({}).as_dict()['detours']
        self.assertEqual( len(def_set), 5)
        sc_scale = 3. * np.identity(3)
        def_task_list = []
        for d in def_set:
            setuptask = d['spec']['_tasks'][0]
            if 'bulk' in d['name']:
                defect_flag = False
                charge = None
            else:
                defect_flag = True
                charge = int( d['name'].split('_')[-2])
                def_task_list.append( setuptask['transformation_params'][0]['defect'])
            self._verify_setuptask( setuptask, supercell = sc_scale,
                                    user_incar_settings={}, charge=charge,
                                    is_defect=defect_flag)

        expected_defects = [ ['Vacancy', None, 'Cs', -1],
                             ['Vacancy', None, 'Cl', 1],
                             ['Substitution', 'Cs', 'Cl', 2],
                             ['Substitution', 'Cl', 'Cs', 0]]
        self._verify_defect_list( def_task_list, expected_defects)


        #check larger supercell, conventional =True, and user_incar_settings changes
        set_incars = {"EDIFF":.000001, "EDIFFG": 0.0001, "ISPIN":1, "LWAVE": False}
        ft = DefectSetupFiretask( structure=struct, cellmax=300, conventional=True,
                                  user_incar_settings= set_incars)

        def_set = ft.run_task({}).as_dict()['detours']
        self.assertEqual( len(def_set), 5)
        sc_scale = 5. * np.identity(3)
        def_task_list = []
        for d in def_set:
            setuptask = d['spec']['_tasks'][0]
            if 'bulk' in d['name']:
                defect_flag = False
                charge = None
            else:
                defect_flag = True
                charge = int( d['name'].split('_')[-2])
                def_task_list.append( setuptask['transformation_params'][0]['defect'])
            self._verify_setuptask( setuptask, supercell = sc_scale,
                                    user_incar_settings=set_incars, charge=charge,
                                    is_defect=defect_flag)

        expected_defects = [ ['Vacancy', None, 'Cs', -1],
                             ['Vacancy', None, 'Cl', 1],
                             ['Substitution', 'Cs', 'Cl', 2],
                             ['Substitution', 'Cl', 'Cs', 0]]
        self._verify_defect_list( def_task_list, expected_defects)

        # Now check different types of defect setups and different types of charges...

        # subcase 1 = default vacancies with custom vacancy charges,
        # default substitutions with custom substitution charges,
        # interstitials with Voronoi site finder and custom charges
        v = []
        subs = {}
        inter = []
        in_chg = {'vacancies': {'Cs': [-3], 'Cl': [2, 1]},
                  'substitutions': {'Cl': {'Cs': [-1]},
                                    'Cs': {'Cl': [4, 3]}}}
        ft = DefectSetupFiretask( structure=struct, cellmax=100,
                                  vacancies=v, substitutions=subs,
                                  interstitials=inter, initial_charges=in_chg)

        def_set = ft.run_task({}).as_dict()['detours']
        self.assertEqual( len(def_set), 7)
        def_task_list = []
        sc_scale = 3. * np.identity(3)
        for d in def_set:
            setuptask = d['spec']['_tasks'][0]
            if 'bulk' in d['name']:
                defect_flag = False
                charge = None
            else:
                defect_flag = True
                charge = int( d['name'].split('_')[-2])
                def_task_list.append( setuptask['transformation_params'][0]['defect'])
            self._verify_setuptask( setuptask, supercell = sc_scale,
                                    user_incar_settings={}, charge=charge,
                                    is_defect=defect_flag)

        expected_defects = [ ['Vacancy', None, 'Cs', -3],
                             ['Vacancy', None, 'Cl', 2],
                             ['Vacancy', None, 'Cl', 1],
                             ['Substitution', 'Cs', 'Cl', -1],
                             ['Substitution', 'Cl', 'Cs', 4],
                             ['Substitution', 'Cl', 'Cs', 3]]
        self._verify_defect_list( def_task_list, expected_defects)

        # subcase 2 = custom vacancies with default charges,
        # custom substitutions (antisite and sub based) with default charges,
        # manual insertion of interstitial defects to use
        v = ['Cl']
        subs = {'Cl': ['Cs'], 'Cs': ['F']}
        defect_site_1 = PeriodicSite( 'H', [0., 1.05225, 2.1045],
                                      struct.lattice,
                                      coords_are_cartesian=True)
        defect_site_2 = PeriodicSite( 'O', [0., 1.05225, 2.1045],
                                      struct.lattice,
                                      coords_are_cartesian=True)
        inter = [['H', Interstitial(struct, defect_site_1, charge=0.)],
                 ['O', Interstitial(struct, defect_site_2, charge=0.)]]
        in_chg = {}
        ft = DefectSetupFiretask( structure=struct, cellmax=100,
                                  vacancies=v, substitutions=subs,
                                  interstitials=inter, initial_charges=in_chg)

        def_set = ft.run_task({}).as_dict()['detours']
        self.assertEqual( len(def_set), 6)
        def_task_list = []
        sc_scale = 3. * np.identity(3)
        for d in def_set:
            setuptask = d['spec']['_tasks'][0]
            if 'bulk' in d['name']:
                defect_flag = False
                charge = None
            else:
                defect_flag = True
                charge = int( d['name'].split('_')[-2])
                def_task_list.append( setuptask['transformation_params'][0]['defect'])
            self._verify_setuptask( setuptask, supercell = sc_scale,
                                    user_incar_settings={}, charge=charge,
                                    is_defect=defect_flag)

        expected_defects = [ ['Vacancy', None, 'Cl', 1],
                             ['Substitution', 'Cs', 'Cl', 2],
                             ['Substitution', 'F', 'Cs', -2],
                             ['Interstitial', 'H', None, 0],
                             ['Interstitial', 'O', None, 0]]
        self._verify_defect_list( def_task_list, expected_defects)

        # subcase 3 = custom vacancies with custom vacancy charges,
        # custom substitutions (antisite and sub based) with custom vacancy charges
        # in custom charges, give an excessive set of charge suggestions for defects
        # that were not requested... this should NOT generate actual defects
        # because we did not ask for those defect types in the v and subs lists
        v = ['Cl']
        subs = {'Cs': ['Cl', 'F']}
        inter = []
        in_chg = {'vacancies': {'Cs': [-3], 'Cl': [3, 4]},
                  'substitutions': {'Cs': {'Cl': [0, -1],
                                           'F': [0]}}}
        ft = DefectSetupFiretask( structure=struct, cellmax=100,
                                  vacancies=v, substitutions=subs,
                                  interstitials=inter, initial_charges=in_chg)

        def_set = ft.run_task({}).as_dict()['detours']
        self.assertEqual( len(def_set), 6)
        def_task_list = []
        sc_scale = 3. * np.identity(3)
        for d in def_set:
            setuptask = d['spec']['_tasks'][0]
            if 'bulk' in d['name']:
                defect_flag = False
                charge = None
            else:
                defect_flag = True
                charge = int( d['name'].split('_')[-2])
                def_task_list.append( setuptask['transformation_params'][0]['defect'])
            self._verify_setuptask( setuptask, supercell = sc_scale,
                                    user_incar_settings={}, charge=charge,
                                    is_defect=defect_flag)

        expected_defects = [ ['Vacancy', None, 'Cl', 3],
                             ['Vacancy', None, 'Cl', 4],
                             ['Substitution', 'Cl', 'Cs', 0],
                             ['Substitution', 'Cl', 'Cs', -1],
                             ['Substitution', 'F', 'Cs', 0]]
        self._verify_defect_list( def_task_list, expected_defects)


        # subcase 4 = custom interstitials with InFit
        # interstitial finder (time intensive) while turning off
        # substitutions and vacancies. Make use of smaller structure to speed up test
        struct = PymatgenTest.get_structure("He_BCC")
        ft = DefectSetupFiretask( structure=struct, cellmax=100,
                                 vacancies=[[]], substitutions={None:{}},
                                 interstitials=[['H', 'InFit']], initial_charges={})

        def_set = ft.run_task({}).as_dict()['detours']
        self.assertEqual( len(def_set), 2)
        def_task_list = []
        for d in def_set:
            setuptask = d['spec']['_tasks'][0]
            if 'bulk' in d['name']:
                defect_flag = False
                charge = None
            else:
                defect_flag = True
                charge = int( d['name'].split('_')[-2])
                def_task_list.append( setuptask['transformation_params'][0]['defect'])
            self._verify_setuptask( setuptask, supercell = sc_scale,
                                    user_incar_settings={}, charge=charge,
                                    is_defect=defect_flag)

        expected_defects = [ ['Interstitial', 'H', None, 0]]
        self._verify_defect_list( def_task_list, expected_defects)


if __name__ == '__main__':
    unittest.main()
