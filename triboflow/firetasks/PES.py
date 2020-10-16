#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 14:52:57 2020

@author: wolloch
"""
import monty
from operator import itemgetter
from pymatgen.core.structure import Structure
from pymatgen.core.operations import SymmOp
from fireworks import FWAction, FiretaskBase, Workflow
from fireworks.utilities.fw_utilities import explicit_serialize
from atomate.utils.utils import env_chk
from atomate.vasp.fireworks.core import OptimizeFW, ScanOptimizeFW
from atomate.vasp.powerups import add_modify_incar
from triboflow.PES_functions.HS_functions import GetSlabHS, GetInterfaceHS
from triboflow.PES_functions.utility_functions import ApplyPbcToHS
from triboflow.helper_functions import GetInterfaceFromDB, \
    GetCustomVaspRelaxSettings, GetDB, GetHighLevelDB, CleanUpSitePorperties

@explicit_serialize
class FT_RetrievePESEnergies(FiretaskBase):
    required_params = ['interface_name', 'functional', 'tag']
    optional_params = ['db_file']
    def run_task(self, fw_spec):
        name = self.get('interface_name')
        functional = self.get('functional')
        tag = self.get('tag')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        
        interface_dict = GetInterfaceFromDB(name, db_file, functional)
        lateral_shifts = interface_dict['PES']['high_symmetry_points']['combined_unique']
        
        DB = GetDB(db_file)
        energy_list=[]
        for s in lateral_shifts.keys():
            label = tag + '_' + s
            x_shift = lateral_shifts.get(s)[0][0]
            y_shift = lateral_shifts.get(s)[0][1]
            vasp_calc = DB.tasks.find_one({'task_label': label})
            energy = vasp_calc['output']['energy']
            energy_list.append([s, x_shift, y_shift, energy])
            
        sorted_energy_list = sorted(energy_list, key=itemgetter(3))
        
        min_stacking = sorted_energy_list[0][0]
        max_stacking = sorted_energy_list[-1][0]
        calc_min = DB.tasks.find_one({'task_label': tag+'_'+min_stacking})
        calc_max = DB.tasks.find_one({'task_label': tag+'_'+max_stacking})
        struct_min = calc_min['output']['structure']
        struct_max = calc_max['output']['structure']
        
        tribo_db = GetHighLevelDB(db_file)
        coll = tribo_db[functional+'.interface_data']
        coll.update_one({'name': name},
                        {'$set': {'relaxed_structrue@min': struct_min,
                                  'relaxed_structrue@max': struct_max,
                                  'PES.high_symmetry_points.energy_list': 
                                         sorted_energy_list}})
        

@explicit_serialize
class FT_FindHighSymmPoints(FiretaskBase):
    required_params = ['interface_name', 'functional']
    optional_params = ['db_file', 'top_name', 'bottom_name']
    def run_task(self, fw_spec):
        name = self.get('interface_name')
        functional = self.get('functional')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        top_name = self.get('top_name', 'top_aligned')
        bottom_name = self.get('bottom_name', 'bottom_aligned')
        
        interface_dict = GetInterfaceFromDB(name, db_file, functional)           
                
        top_aligned = Structure.from_dict(interface_dict[top_name])
        bottom_aligned = Structure.from_dict(interface_dict[bottom_name])
        
        #top slab needs to be mirrored to find the high symmetry points at the
        #interface.
        mirror = SymmOp.reflection(normal=[0,0,1], origin=[0, 0, 0])
        flipped_top = top_aligned.copy()
        flipped_top.apply_operation(mirror)
        top_hsp_unique, top_hsp_all = GetSlabHS(flipped_top)
        
        bottom_hsp_unique, bottom_hsp_all = GetSlabHS(bottom_aligned)
        
        hsp_unique = GetInterfaceHS(bottom_hsp_unique, top_hsp_unique)
        hsp_all = GetInterfaceHS(bottom_hsp_all, top_hsp_all)
        
# =============================================================================
#       Project the high symmetry points which might be outside of the cell
#       back into it. SHOULD DEFINATELY BE REVISED WITHOUT USING THE PYTORCH
#       PACKAGE!
# =============================================================================
        hsp_unique = ApplyPbcToHS(bottom_aligned, hsp_unique)
        hsp_all = ApplyPbcToHS(bottom_aligned, hsp_all)
        
        b_hsp_u =  monty.json.jsanitize(bottom_hsp_unique)
        b_hsp_a =  monty.json.jsanitize(bottom_hsp_all)
        t_hsp_u =  monty.json.jsanitize(top_hsp_unique)        
        t_hsp_a =  monty.json.jsanitize(top_hsp_all)
        c_hsp_u =  monty.json.jsanitize(hsp_unique)        
        c_hsp_a =  monty.json.jsanitize(hsp_all)
        
        tribo_db = GetHighLevelDB(db_file)
        coll = tribo_db[functional+'.interface_data']
        
        coll.update_one({'name': name},
                        {'$set': {'PES':
                                    {'high_symmetry_points': 
                                         {'bottom_unique': b_hsp_u,
                                          'bottom_all': b_hsp_a,
                                          'top_unique': t_hsp_u,
                                          'top_all': t_hsp_a,
                                          'combined_unique': c_hsp_u,
                                          'combined_all': c_hsp_a}}}})
            
        return FWAction(update_spec=({'lateral_shifts': c_hsp_u}))
        
        

@explicit_serialize
class FT_StartPESCalcs(FiretaskBase):
    """Start z-relaxations for different lateral positions of an interface.
    
    Take a list of lateral shifts from the fw_spec and start relaxations
    for each one of them as parallel detours.
    
    Parameters
    ----------
    interface_name : str
        Name of the interface in the high-level database.
    functional : str
        Which functional to use; has to be 'PBE' or 'SCAN'.
    tag : str
        Unique tag to identify the calculations.
    db_file : str
        Full path to the db.json file that should be used. Defaults to
        '>>db_file<<', to use env_chk.
    structure_name : str, optional
        Name of the structure in the interface entry to the high-level database
        for which the PPES should be calculated. The default is
        'unrelaxed_structure'.
        
    Returns
    -------
    FWActions that produce a detour workflow with relaxations for the PES.
    """
    required_params = ['interface_name', 'functional', 'tag']
    optional_params = ['db_file', 'structure_name']
    def run_task(self, fw_spec):
        name = self.get('interface_name')
        functional = self.get('functional')
        tag = self.get('tag')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        structure_name = self.get('structure_name', 'unrelaxed_structure')
        
        lateral_shifts = fw_spec.get('lateral_shifts')
        if not lateral_shifts:
            raise SystemExit('Lateral shifts not found in the fw_spec./n'
                             'Please check your Firework for errors!')
        
        interface_dict = GetInterfaceFromDB(name, db_file, functional)            
                
        comp_params = interface_dict['comp_parameters']
        struct = Structure.from_dict(interface_dict[structure_name])
        
        # List all sites of interface that have positive c coordinates as they
        # are in the upper slab.
        sites_to_shift = []
        for i, s in enumerate(struct.sites):
            if s.c > 0:
                sites_to_shift.append(i)
        
        FW_list=[]
        for s in lateral_shifts.keys():
            label = tag + '_' + s
            x_shift = lateral_shifts.get(s)[0][0]
            y_shift = lateral_shifts.get(s)[0][1]
            #Make sure that there are no NoneTypes in the site_properties!
            struct_s = CleanUpSitePorperties(struct.copy())
            struct_s.translate_sites(indices=sites_to_shift,
                                     vector=[x_shift, y_shift, 0],
                                     frac_coords=False, to_unit_cell=False)
            
            vis = GetCustomVaspRelaxSettings(structure=struct_s,
                                             comp_parameters=comp_params,
                                             relax_type='interface_z_relax')
            if functional == 'SCAN':
                FW = ScanOptimizeFW(structure=struct_s,
                                    name=label,
                                    vasp_input_set = vis)
            else:
                FW = OptimizeFW(structure=struct_s,
                                name=label,
                                vasp_input_set = vis)
            FW_list.append(FW)
            
        
        WF = Workflow(FW_list, name='PES relaxations for: '+name)
        PES_Calcs_WF = add_modify_incar(WF)
        
        return FWAction(detours = PES_Calcs_WF)