#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 16:11:53 2021

@author: mwo
"""

from fireworks import FWAction, FiretaskBase, Firework, Workflow
from fireworks.utilities.fw_utilities import explicit_serialize

from atomate.vasp.fireworks.core import StaticFW

from triboflow.utils.database import StructureNavigator, Navigator
from triboflow.utils.vasp_tools import get_custom_vasp_static_settings

@explicit_serialize
class FT_GetEpsilon(FiretaskBase):
    _fw_name = 'Read the static dielectric constant'
    required_params = ['label']
    optional_params = ['db_file']
    def run_task(self, fw_spec):

        label = self.get('label')
        db_file = self.get('db_file', 'auto')
        
        
        nav = Navigator(db_file)
        output_data = nav.find_data(collection='tasks',
                                    fltr={'task_label': label})
        eps_tensor = output_data['output']['epsilon_static']
        eps_average = sum([eps_tensor[0][0], eps_tensor[1][1], eps_tensor[1][1]])/3.0
        
        return FWAction(update_spec={'epsilon_tensor': eps_tensor,
                                     'epsilon': eps_average})

@explicit_serialize
class FT_UpdateCompParams(FiretaskBase):
    _fw_name = 'Update computational parameters'
    required_params = ['mpid', 'functional', 'new_params']
    optional_params = ['db_file', 'update_bulk', 'update_slabs',
                       'high_level_db']
    def run_task(self, fw_spec):

        mpid = self.get('mpid')
        functional = self.get('functional')
        new_params = self.get('new_params')
        update_bulk = self.get('update_bulk', True)
        update_slabs = self.get('update_bulk', False)
        db_file = self.get('db_file', 'auto')
        hl_db = self.get('high_level_db', True)
        
        nav_structure = StructureNavigator(db_file=db_file, 
                                           high_level=hl_db)
        
        #get values from spec:
        new_data = {}
        for param in new_params:
            new_data[param] = fw_spec.get(param, None)
            
        
        if update_bulk:
            nav_structure.update_data(collection=functional+'.bulk_data',
                                      fltr={'mpid': mpid},
                                      new_values=new_data,
                                      upsert=False)
        if update_slabs:
            nav_structure.update_data(collection=functional+'.slab_data',
                                      fltr={'mpid': mpid},
                                      new_values=new_data,
                                      upsert=False)
        return FWAction(update_spec=fw_spec)
        

def dielectric_constant_swf(structure,
                            mpid,
                            flag, 
                            comp_parameters={}, 
                            spec={},
                            functional='PBE', 
                            db_file=None,
                            hl_db=True,
                            update_bulk=True,
                            update_slabs=False):
    
    formula = structure.composition.reduced_formula
    wf_name = f'Dielectric constant calculation workflow for {formula}'
    
    vis = get_custom_vasp_static_settings(structure,
                                          comp_parameters,
                                          'bulk_epsilon_from_scratch',
                                          name=flag)
    Calc_Eps_FW = StaticFW(structure=structure,
                           name=flag,
                           vasp_input_set=vis)
    
    Get_Eps_FT = FT_GetEpsilon(label=flag, db_file=db_file)
    Update_Data_FT = FT_UpdateCompParams(mpid=mpid,
                                         functional=functional,
                                         new_params='epsilon',
                                         update_bulk=update_bulk,
                                         update_slabs=update_slabs,
                                         db_file=db_file,
                                         hl_db=hl_db)
    Update_FW = Firework(tasks=[Get_Eps_FT, Update_Data_FT],
                         spec=spec,
                         name=flag+'_update_high_level')
    
    WF = Workflow([Calc_Eps_FW, Update_FW], {Calc_Eps_FW: [Update_FW]},
                  name=wf_name)
    return WF