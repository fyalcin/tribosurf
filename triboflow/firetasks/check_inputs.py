"""Utility Firetasks for the triboflow project.

Created on Mon Jun 22 12:29:28 2020
@author: mwo
"""
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from fireworks import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize
from atomate.utils.utils import env_chk

from triboflow.utils.database import Navigator, StructureNavigator, NavigatorMP
from triboflow.utils.structure_manipulation import interface_name


@explicit_serialize
class FT_UpdateCompParams(FiretaskBase):

    _fw_name = 'Update computational parameters in high level DB'
    required_params = ['mp_id_1', 'miller_1', 'mp_id_2', 'miller_2',
                       'functional']
    optional_params = ['db_file']

    def run_task(self, fw_spec):

        mp_id_1 = self.get('mp_id_1')
        mp_id_2 = self.get('mp_id_2')

        if type(self['miller_1']) == str:
            miller_1 = [int(k) for k in list(self['miller_1'])]
        else:
            miller_1 = self['miller_1']
        if type(self['miller_2']) == str:
            miller_2 = [int(k) for k in list(self['miller_2'])]
        else:
            miller_2 = self['miller_2']
        functional = self.get('functional')

        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        
        nav_structure = StructureNavigator(
            db_file=db_file, 
            high_level='triboflow')
        bulk_1 = nav_structure.get_bulk_from_db(
            mp_id=mp_id_1, 
            functional=functional)
        bulk_2 = nav_structure.get_bulk_from_db(
            mp_id=mp_id_2,
            functional=functional
        )

        encut_1 = bulk_1['comp_parameters']['encut']
        encut_2 = bulk_2['comp_parameters']['encut']
        encut_inter = max(encut_1, encut_2)
        k_dens_1 = bulk_1['comp_parameters']['k_dens']
        k_dens_2 = bulk_2['comp_parameters']['k_dens']
        k_dens_inter = max(k_dens_1, k_dens_2)
        metal_1 = bulk_1['comp_parameters']['is_metal']
        metal_2 = bulk_2['comp_parameters']['is_metal']
        metal_inter = any((metal_1, metal_2))

        nav_high = Navigator(db_file=db_file, high_level='triboflow')
        nav_high.update_data(
            collection=functional+'.slab_data',
            fltr={'mpid': mp_id_1, 'miller': miller_1},
            new_values={'$set': {'comp_parameters.encut': encut_1,
                                 'comp_parameters.k_dens': k_dens_1}})
        nav_high.update_data(
            collection=functional+'.slab_data',
            fltr={'mpid': mp_id_2, 'miller': miller_2},
            new_values={'$set': {'comp_parameters.encut': encut_2,
                                 'comp_parameters.k_dens': k_dens_2}})

        inter_name = interface_name(mp_id_1, miller_1, mp_id_2, miller_2)
        nav_high.update_data(
            collection=functional+'.interface_data',
            fltr={'name': inter_name},
            new_values={'$set': {'comp_parameters.encut': encut_inter,
                                 'comp_parameters.k_dens': k_dens_inter,
                                 'comp_parameters.is_metal': metal_inter}})

@explicit_serialize
class FT_MakeInterfaceInDB(FiretaskBase):

    _fw_name = 'Make bulk entry into high level DB'
    required_params = ['mat1_data_loc', 'mat2_data_loc', 'comp_data_loc',
                       'interface_data_loc']
    optional_params = ['db_file']

    def run_task(self, fw_spec):
        data1 = fw_spec[self['mat1_data_loc']]
        data2 = fw_spec[self['mat2_data_loc']]
        comp_data = fw_spec[self['comp_data_loc']]
        interface_data = fw_spec[self['interface_data_loc']]

        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        
        functional = comp_data['functional']

        nav_mp = NavigatorMP()
        struct1, mp_id_1 = nav_mp.get_low_energy_structure(
            chem_formula=data1['formula'],
            mp_id=data1['mp_id'])
        struct2, mp_id_2 = nav_mp.get_low_energy_structure(
            chem_formula=data2['formula'],
            mp_id=data2['mp_id'])

        nav_high = Navigator(db_file=db_file, high_level='triboflow')

        name = interface_name(mp_id_1, data1['miller'],
                              mp_id_2, data2['miller'])

        if nav_high.find_data(
               collection=functional+'.interface_data', 
               fltr={'name': name}):

            print('{} interface can not be added to {}.interface_data '
                  'collection because an interface with that name is already '
                  'present.'.format(name, functional))
            return
        else:
            nav_high.insert_data(
                collection=functional+'.interface_data',
                data={'name': name,
                      'comp_parameters': comp_data,
                      'interface_parameters': interface_data})
            return

@explicit_serialize
class FT_MakeSlabInDB(FiretaskBase):

    _fw_name = 'Make bulk entry into high level DB'
    required_params = ['mat_data_loc', 'comp_data_loc']
    optional_params = ['db_file']

    def run_task(self, fw_spec):

        data = fw_spec[self['mat_data_loc']]
        comp_data = fw_spec[self['comp_data_loc']]

        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        
        functional = comp_data['functional']
        
        nav_mp = NavigatorMP()
        struct, mp_id = nav_mp.get_low_energy_structure(
            chem_formula=data['formula'], 
            mp_id=data['mp_id'])

        bandgap = nav_mp.get_property_from_mp(
            mp_id=mp_id, 
            properties=['band_gap'])
        bandgap = bandgap['band_gap']

        if bandgap > 0.3:
            comp_data['is_metal'] = False
        else:
            comp_data['is_metal'] = True

        nav_high = Navigator(db_file=db_file, high_level='triboflow')

        if nav_high.find_data(
            collection=functional+'.slab_data',
            fltr={'mpid': mp_id, 'miller': data['miller']}):

            print('{}-{} slab can not be added to {}.slab_data collection'
                  'because a material with MP-ID {} is already present.'
                  .format(data['formula'], data['miller'], mp_id, functional))
            return
        else:
            nav_high.insert_data(
                collection=functional+'.slab_data',
                data={'mpid': mp_id,
                      'formula': data['formula'],
                      'miller': data['miller'],
                      'min_thickness': data['min_thickness'],
                      'min_vacuum': data['min_vacuum'],
                      'comp_parameters': comp_data})
            return


@explicit_serialize
class FT_MakeBulkInDB(FiretaskBase):

    _fw_name = 'Make bulk entry into high level DB'
    required_params = ['mat_data_loc', 'comp_data_loc']
    optional_params = ['db_file']

    def run_task(self, fw_spec):

        data = fw_spec[self['mat_data_loc']]
        comp_data = fw_spec[self['comp_data_loc']]

        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        
        functional = comp_data['functional']

        nav_mp = NavigatorMP()
        struct, mp_id = nav_mp.get_low_energy_structure(
            chem_formula=data['formula'], 
            mp_id=data['mp_id'])
        sga = SpacegroupAnalyzer(struct)
        prim_struct = sga.get_primitive_standard_structure()

        bandgap = nav_mp.get_property_from_mp(
            mp_id=mp_id, 
            properties=['band_gap'])
        bandgap = bandgap['band_gap']

        if bandgap > 0.3:
            comp_data['is_metal'] = False
        else:
            comp_data['is_metal'] = True
        
        nav_high = Navigator(db_file=db_file, high_level='triboflow')

        if nav_high.find_data(
            collection=functional+'.bulk_data',
            fltr={'mpid': mp_id}):

            print('{} bulk can not be added to bulk_data collection because a'
                  'material with MP-ID {} is already present in the {} '
                  'collection!'.format(data['formula'], mp_id, functional))
            return
        else:
            nav_high.insert_data(
                collection=functional+'.bulk_data',
                data={'mpid': mp_id,
                        'formula': data['formula'],
                        'structure_fromMP': struct.as_dict(),
                        'primitive_structure': prim_struct.as_dict(),
                        'comp_parameters': comp_data})
            return
        
    

@explicit_serialize
class FT_CheckCompParamDict(FiretaskBase):
    """Checks a dictionary for essential keys and adds default values if needed.
    
    An input dictionary is compared to a list of essential keys that must be
    present in the dictionary. If essential keys are missing an error is
    printed and SystemExit is raised. If optional keys are not given, default
    values are used for them and added to the dictionary.
    
    Parameters
    ----------
    input_dict_loc: list of str
        The location of a dictionary in the spec that holds the computational
        input parameters
    output_dict_loc: list of str, optional
        The location of the output dictionary that is going to be put into the
        fw_spec. If not specified this will default to the input_dict_loc.
    
    Returns
    -------
    dict
        A dictionary is pushed to the fw_spec. It contains all the
        essential and optional keys given and uses default values for the
        remaining optional keys.
    """

    _fw_name = 'Check Input Dict'
    required_params = ['input_dict']
    optional_params = ['output_dict_name']
    
    def run_task(self, fw_spec):
        """Run the FireTask.

        Edit the essential and additional keys here if necessary, as well as
        the default values given for the additional keys.
        """
        #Edit this block according to need, but be careful with the defaults!
        #####################################################################
        essential_keys = ['use_vdw']
        additional_keys = ['volume_tolerance', 'energy_tolerance', 'use_spin',
                           'BM_tolerance', 'functional']
        
        volume_tolerance_default = 0.001
        energy_tolerance_default = 0.001
        functional_default = 'SCAN'
        use_spin_default = True
        BM_tolerance_default = 0.01
        #####################################################################
        
        input_dict = self['input_dict']
        output_dict_name = self.get('output_dict_name', 'comp_parameters')

        # Define all known keys here
        known_keys = essential_keys + additional_keys

        # Create lists for input options considered to be True:
        true_list = ['true', 'True', 'TRUE', '.TRUE.', '.true.', True, 
                     'yes', 'Yes', 'YES', '.YES.', '.yes.']

        # Initialize checking dictionary for essential inputs and output dict
        out_dict = {}
        check_essential = {}
        for key in essential_keys:
            check_essential[key] = False
        for key in input_dict.keys():
        # Check for unknown parameters:
            if key not in known_keys:
                raise SystemExit('The input parameter <'+str(key)+
                              '> is not known. Please check your input file'
                              'and use only the following parameters:\n'
                              '{}'.format(known_keys))
            elif key == 'use_vdw':
                if input_dict[key] in true_list:
                    out_dict['use_vdw'] = True
                else:
                    out_dict['use_vdw'] = False
                check_essential[key] = True
        
        
        # Check if all essential keys are present and print out missing ones.
        missing_keys=[]
        for key in check_essential.keys():
            if check_essential[key] == False:
                missing_keys.append(key)
        if not all(check_essential.values()):
            raise SystemExit('The essential input parameters "'+missing_keys+
                              '" are missing. Check your input file!\n')
        
        for key in additional_keys:
            if key == 'use_spin':
                if key in input_dict:
                    if input_dict[key] in true_list:
                        out_dict['use_spin'] = True
                    else:
                        out_dict['use_spin'] = False
                else:
                    out_dict['use_spin'] = use_spin_default
            if key == 'volume_tolerance':
                if key in input_dict:
                    out_dict['volume_tolerance'] = float(input_dict[key])
                else:
                    out_dict['volume_tolerance'] = volume_tolerance_default
            if key == 'energy_tolerance':
                if key in input_dict:
                    out_dict['energy_tolerance'] = float(input_dict[key])
                else:
                    out_dict['energy_tolerance'] = energy_tolerance_default
            if key == 'functional':
                if key in input_dict:
                    out_dict['functional'] = str(input_dict[key])
                else:
                    out_dict['functional'] = functional_default
            if key == 'BM_tolerance':
                if key in input_dict:
                    out_dict['BM_tolerance'] = float(input_dict[key])
                else:
                    out_dict['BM_tolerance'] = BM_tolerance_default
                
        return FWAction(mod_spec=[{'_set': {output_dict_name: out_dict}}])
    
@explicit_serialize
class FT_CheckInterfaceParamDict(FiretaskBase):
    """Checks a dictionary for essential keys and adds default values if needed.
    
    An input dictionary is compared to a list of essential keys that must be
    present in the dictionary. If essential keys are missing an error is
    printed and SystemExit is raised. If optional keys are not given, default
    values are used for them and added to the dictionary.
    
    Parameters
    ----------
    input_dict_loc: list of str
        The location of a dictionary in the spec that holds the computational
        input parameters
    output_dict_loc: list of str, optional
        The location of the output dictionary that is going to be put into the
        fw_spec. If not specified this will default to the input_dict_loc.
    
    Returns
    -------
    dict
        A dictionary is pushed to the fw_spec. It contains all the
        essential and optional keys given and uses default values for the
        remaining optional keys.
    """

    _fw_name = 'Check Input Dict'
    required_params = ['input_dict']
    optional_params = ['output_dict_name']
    
    def run_task(self, fw_spec):
        """Run the FireTask.

        Edit the essential and additional keys here if necessary, as well as
        the default values given for the additional keys.
        """
        #Edit this block according to need, but be careful with the defaults!
        #####################################################################
        essential_keys = ['max_area']
        additional_keys = ['interface_distance', 'max_mismatch',
                           'max_angle_diff', 'r1r2_tol']
        
        defaults = {'interface_distance': 2.0,
                    'max_mismatch': 0.01,
                    'max_angle_diff': 1.5,
                    'r1r2_tol': 0.05}
        #####################################################################
         
        input_dict = self['input_dict']
        output_dict_name = self.get('output_dict_name', 'interface_params')

        # Define all known keys here
        known_keys = essential_keys + additional_keys

        # Initialize checking dictionary for essential inputs and output dict
        out_dict = {}
        check_essential = {}
        for key in essential_keys:
            check_essential[key] = False
        for key in input_dict.keys():
        # Check for unknown parameters:
            if key not in known_keys:
                raise SystemExit('The input parameter <'+str(key)+
                              '> is not known. Please check your input file'
                              'and use only the following parameters:\n{}'
                              .format(known_keys))
            elif key == 'max_area':
                out_dict['max_area'] = float(input_dict[key])
                check_essential[key] = True
        
        
        # Check if all essential keys are present and print out missing ones.
        missing_keys=[]
        for key in check_essential.keys():
            if check_essential[key] == False:
                missing_keys.append(key)
        if not all(check_essential.values()):
            raise SystemExit('The essential input parameters "'+missing_keys+
                              '" are missing. Check your input file!\n')
        
        for key in additional_keys:
            if key in input_dict:
                out_dict[key] = float(input_dict[key])
            else:
                out_dict[key] = defaults[key]
                
        return FWAction(mod_spec=[{'_set': {output_dict_name: out_dict}}])

@explicit_serialize
class FT_CheckMaterialInputDict(FiretaskBase):
    """Checks a dictionary for essential keys and adds default values if needed.
    
    An input dictionary is compared to a list of essential keys that must be
    present in the dictionary. If essential keys are missing an error is
    printed and SystemExit is raised. If optional keys are not given, default
    values are used for them and added to the dictionary.
    
    Parameters
    ----------
    input_dict_name: list of str
        List of keys that specify the location of the input parameters
        dictionary in the spec
    output_dict_name: list of str, optional
        Location of the output dictionary that is going to be put into the
        spec. If not specified this will default to the input_dict_loc.
    
    Returns
    -------
    dict
        A dictionary is pushed to the spec of the workflow. It contains all the
        essential and optional keys given and uses default values for the
        remaining optional keys. The name in the spec is output_dict_name.
    """

    _fw_name = 'Check Material Input Dict'
    required_params = ['input_dict']
    optional_params = ['output_dict_name']
    
    def run_task(self, fw_spec):
        """Run the FireTask."""
        #Edit this block according to need, but be careful with the defaults!
        #####################################################################
        essential_keys = ['formula', 'miller']
        additional_keys = ['min_vacuum', 'mp_id', 'min_thickness']
        
        min_thickness_default = 10.0
        min_vacuum_default = 25.0
        # MPID of the minimum energy structure for this formula will be used
        # as default.
        #####################################################################
        
        input_dict = self['input_dict']
        output_dict_name = self.get('output_dict_name', 'materials_params')

        # Define all known keys here
        known_keys = essential_keys + additional_keys

        # Initialize checking dictionary for essential inputs and output dict
        out_dict = {}
        check_essential = {}
        for key in essential_keys:
            check_essential[key] = False
            
        for key in input_dict.keys():
        # Check for unknown parameters:
            if key not in known_keys:
                with open('output.txt', 'a') as out:
                    out.write(' ')
                    out.write('The input parameter <'+str(key)+
                              '> is not known. Please check your input file'
                              'and use only the following parameters:\n')
                    out.write(str(known_keys))
                    out.write('\n')
                raise SystemExit
            elif key == 'formula':
                out_dict['formula'] = str(input_dict[key])
                check_essential[key] = True
            elif key == 'miller':
                out_dict['miller'] = [int(k) for k in list(input_dict[key])]
                check_essential[key] = True
        
        # Check if all essential keys are present and print out missing ones.
        for key in check_essential.keys():
            if check_essential[key] == False:
                with open('output.txt', 'a') as out:
                    out.write('The essential input parameter "'+key+
                              '" is missing. Check your input file!\n')
        if not all(check_essential.values()):
            with open('output.txt', 'a') as out:
                out.write('')
            raise SystemExit
            
        for key in additional_keys:
            if key == 'mp_id':
                if key in input_dict:
                    out_dict['mp_id'] = str(input_dict[key])
                else:
                    nav_mp = NavigatorMP()
                    struct, mp_id = nav_mp.get_low_energy_structure(
                        chem_formula=str(input_dict['formula']))
                    out_dict['mp_id'] = mp_id
            
            if key == 'min_thickness':
                if key in input_dict:
                    out_dict['min_thickness'] = float(input_dict[key])
                else:
                    out_dict['min_thickness'] = min_thickness_default
            
            if key == 'min_vacuum':
                if key in input_dict:
                    out_dict['min_vacuum'] = float(input_dict[key])
                else:
                    out_dict['min_vacuum'] = min_vacuum_default
        
        return FWAction(mod_spec=[{'_set': {output_dict_name: out_dict}}])

