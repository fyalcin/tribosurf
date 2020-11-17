#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 14:52:57 2020

@author: wolloch
"""
import monty
import numpy as np
from operator import itemgetter
from scipy.interpolate import Rbf
from pymatgen.core.structure import Structure
from pymatgen.core.surface import Slab
from pymatgen.core.operations import SymmOp
from fireworks import FWAction, FiretaskBase, Workflow
from fireworks.utilities.fw_utilities import explicit_serialize
from atomate.utils.utils import env_chk
from atomate.vasp.fireworks.core import OptimizeFW, ScanOptimizeFW
from atomate.vasp.powerups import add_modify_incar
from triboflow.phys.high_symmetry import GetSlabHS, GetInterfaceHS, \
    PBC_HSPoints, RemoveDuplicatesFromHSDicts
from triboflow.phys.potential_energy_surface import GetPES
from triboflow.utils.database import GetInterfaceFromDB, GetDB, GetHighLevelDB
from triboflow.utils.vasp_tools import GetCustomVaspRelaxSettings
from triboflow.utils.structure_manipulation import CleanUpSiteProperties, \
    SlabFromStructure


@explicit_serialize
class PreparePesCompute_FT(FiretaskBase):
    required_params = ['interface_name', 'functional']
    optional_params = ['db_file']
    def run_task(self, fw_spec):
        name = self.get('interface_name')
        functional = self.get('functional')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        
        
        inter_dict = GetInterfaceFromDB(name, db_file, functional)
        struct = Structure.from_dict(inter_dict['relaxed_structure@min'])
        
        #Copy the energies for the unique points to all points
        E_unique = inter_dict['PES']['high_symmetry_points']['energy_list']
        all_hs = inter_dict['PES']['high_symmetry_points']['combined_all']
        
        Interpolation, pes_dict, pes_data = GetPES(hs_all=all_hs,
                                                   E=E_unique,
                                                   cell=struct.lattice.matrix,
                                                   to_fig=False)
        
        
        
        
        

@explicit_serialize
class FT_RetrievePESEnergies(FiretaskBase):
    """Retrieve the energies from the PES relaxations and update the db.
    
    Uses a tag together with the labels of the high-symmetry points saved
    in the high level database to retrieve the correct energies for each
    lateral shift of the interface. Sort the shifts by energies and save both
    the configuration with the lowest and the highest energy in the high level
    database. Also save the list of shifts and energies with corresponding
    labels there.

    Parameters
    ----------
    interface_name : str
        Name of the interface in the high-level database.
    functional : str
        Which functional to use; has to be 'PBE' or 'SCAN'.
    tag : str
        Unique tag to identify the calculations.
    db_file : str, optional
        Full path to the db.json file that should be used. Defaults to
        '>>db_file<<', to use env_chk.
        
    Returns
    -------
    FWActions that produce a detour workflow with relaxations for the PES.
    """
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
                        {'$set': {'relaxed_structure@min': struct_min,
                                  'relaxed_structure@max': struct_max,
                                  'PES.high_symmetry_points.energy_list': 
                                         sorted_energy_list}})
        
        

@explicit_serialize
class FT_FindHighSymmPoints(FiretaskBase):
    """Compute high symmetry points for the top and bottom slab and the interface.
    
    Finds the high symmetry points of the top side of the bottom slab and the
    bottom side of the top slab. This is done twice, once omitting duplicates,
    and once allowing them. It is made sure that the results are cartesian
    coordinates that lie inside the unit cell. The lists are combined so that
    every combination of unique points for the interface is present. The PES
    section of the high level database is updated with the results and the
    fw_spec is updated with the lateral shifts needed for the PES relaxations
    as well.
    
    Parameters
    ----------
    interface_name : str
        Name of the interface in the high-level database.
    functional : str
        Which functional to use; has to be 'PBE' or 'SCAN'.
    db_file : str, optional
        Full path to the db.json file that should be used. Defaults to
        '>>db_file<<', to use env_chk.
    top_name : str, optional
        Name of the structure in the interface entry to the high-level database
        which constitutes the upper slab. The default is 'top_aligned'.
    bottom_name : str, optional
        Name of the structure in the interface entry to the high-level database
        which constitutes the lower slab. The default is 'bottom_aligned'.
        
    Returns
    -------
    FWActions that updates the fw_spec with lateral shifts.
    """
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
                
        top_aligned = Slab.from_dict(interface_dict[top_name])
        bottom_aligned = Slab.from_dict(interface_dict[bottom_name])
        
        #top slab needs to be mirrored to find the high symmetry points at the
        #interface.
        mirror = SymmOp.reflection(normal=[0,0,1], origin=[0, 0, 0])
        flipped_top = top_aligned.copy()
        flipped_top.apply_operation(mirror)
        top_hsp_unique, top_hsp_all = GetSlabHS(flipped_top)
        
        bottom_hsp_unique, bottom_hsp_all = GetSlabHS(bottom_aligned)
        
        cell = bottom_aligned.lattice.matrix
        
        hsp_unique = GetInterfaceHS(bottom_hsp_unique, top_hsp_unique, cell)
        hsp_all = GetInterfaceHS(bottom_hsp_all, top_hsp_all, cell)
        
        c_hsp_u, c_hsp_a = RemoveDuplicatesFromHSDicts(hsp_unique,
                                                       hsp_all,
                                                       decimals=5)
        
        cell = bottom_aligned.lattice.matrix
           
        b_hsp_u =  PBC_HSPoints(bottom_hsp_unique, cell)
        b_hsp_a =  PBC_HSPoints(bottom_hsp_all, cell)
        t_hsp_u =  PBC_HSPoints(top_hsp_unique, cell)
        t_hsp_a =  PBC_HSPoints(top_hsp_all, cell)
        
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
    Heterogeneous_WF
    Parameters
    ----------
    interface_name : str
        Name of the interface in the high-level database.
    functional : str
        Which functional to use; has to be 'PBE' or 'SCAN'.
    tag : str
        Unique tag to identify the calculations.
    db_file : str, optional
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
            struct_s = CleanUpSiteProperties(struct.copy())
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