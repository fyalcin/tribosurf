"""Utility Firetasks for the triboflow project.

Created on Mon Jun 22 12:29:28 2020
@author: mwo
"""

from fireworks import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize
from atomate.utils.utils import env_chk
from triboflow.helper_functions import GetBulkFromD


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

        #Define all known keys here
        known_keys = essential_keys + additional_keys
        #Create lists for input options condidered to be True:
        true_list = ['true', 'True', 'TRUE', '.TRUE.', '.true.', True, 
                     'yes', 'Yes', 'YES', '.YES.', '.yes.']
        #initialize checking dictionary for essential inputs and output dict
        out_dict = {}
        check_essential = {}
        for key in essential_keys:
            check_essential[key] = False
        for key in input_dict.keys():
        #Check for unknown parameters:
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
        
        
        #check if all essential keys are present and print out missing ones.
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
    required_params = ['input_dict_loc']
    optional_params = ['output_dict_loc']
    
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
                    'max_angle_diff': 2.0,
                    'r1r2_tol': 0.05}
        #####################################################################
        
        input_dict = GetValueFromNestedDict(fw_spec, self['input_dict_loc'])
        
        if 'output_dict_loc' in self:
            out_loc = self['output_dict_loc']
        else:
            out_loc = self['input_dict_loc']

        #Define all known keys here
        known_keys = essential_keys + additional_keys
        #initialize checking dictionary for essential inputs and output dict
        out_dict = {}
        check_essential = {}
        for key in essential_keys:
            check_essential[key] = False
        for key in input_dict.keys():
        #Check for unknown parameters:
            if key not in known_keys:
                raise SystemExit('The input parameter <'+str(key)+
                              '> is not known. Please check your input file'
                              'and use only the following parameters:\n'
                              .format(known_keys))
            elif key == 'max_area':
                out_dict['max_area'] = float(input_dict[key])
                check_essential[key] = True
        
        
        #check if all essential keys are present and print out missing ones.
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
                
        spec = fw_spec
        final_output = WriteNestedDictFromList(out_loc, out_dict, d={})
        updated_spec = UpdateNestedDict(spec, final_output)
        return FWAction(update_spec=updated_spec)

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
        List of keys that specifiy the location of the input parameters
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
    required_params = ['input_dict_loc']
    optional_params = ['output_dict_loc']
    
    def run_task(self, fw_spec):
        """Run the FireTask."""
        #Edit this block according to need, but be careful with the defaults!
        #####################################################################
        essential_keys = ['formula', 'miller']
        additional_keys = ['min_vacuum', 'mp_id', 'min_thickness']
        
        MP_ID_default = None
        min_thickness_default = 10.0
        min_vacuum_default = 25.0
        #####################################################################
        input_dict = GetValueFromNestedDict(fw_spec, self['input_dict_loc'])
        if 'output_dict_name' in self:
            out_loc = self['output_dict_loc']
        else:
            out_loc = self['input_dict_loc']

        #Define all known keys here
        known_keys = essential_keys + additional_keys
        #initialize checking dictionary for essential inputs and output dict
        out_dict = {}
        check_essential = {}
        for key in essential_keys:
            check_essential[key] = False
            
        for key in input_dict.keys():
        #Check for unknown parameters:
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
        
        #check if all essential keys are present and print out missing ones.
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
                    out_dict['mp_id'] = MP_ID_default
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
        
        spec = fw_spec
        spec_update = WriteNestedDictFromList(out_loc, out_dict, d={})
        updated_spec = UpdateNestedDict(spec, spec_update)
        return FWAction(update_spec=updated_spec)

