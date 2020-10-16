#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 16:10:28 2020

@author: mwo
"""

import numpy as np
from scipy.optimize import curve_fit
from pymatgen.core.structure import Structure
from fireworks import FWAction, FiretaskBase, Workflow
from fireworks.utilities.fw_utilities import explicit_serialize
from atomate.utils.utils import env_chk
from atomate.vasp.fireworks.core import StaticFW
from atomate.vasp.powerups import add_modify_incar
from triboflow.helper_functions import GetInterfaceFromDB, \
    GetCustomVaspStaticSettings, GetDB, GetHighLevelDB, CleanUpSitePorperties


@explicit_serialize
class FT_DoPPESCalcs(FiretaskBase):
    """
    Start a CalcPPES_SWF subworkflow that calculates a PPES.
    
    The workflow is only added if there are not already relevant results in
    the high-level database.
    
    Parameters
    ----------
    interface_name : str
        Name of the interface in the high-level database.
    functional : str
        Which functional to use; has to be 'PBE' or 'SCAN'.
    distance_list : list of float, optional
        Modification of the equilibrium distance between the slabs.
        The default is [-0.5, -0.25, 0.0, 0.25, 0.5, 2.5, 3.0, 4.0, 5.0, 7.5].
    out_name : str, optional
        Name for the PPES data in the high-level database. The default is
        'PPES@minimum'.
    structure_name : str, optional
        Name of the structure in the interface entry to the high-level database
        for which the PPES should be calculated. The default is
        'minimum_relaxed'.
    spec : dict, optional
        fw_spec that can be passed to the SWF and will be passed on. The
        default is {}.

    Returns
    -------
    SWF : fireworks.core.firework.Workflow
        Subworkflow to calculate the PPES for a certain interface.

    """
    required_params = ['interface_name', 'functional', 'distance_list']
    optional_params = ['db_file', 'structure_name', 'out_name']
    def run_task(self, fw_spec):
        from triboflow.workflows.subworkflows import CalcPPES_SWF
        name = self.get('interface_name')
        functional = self.get('functional')
        tag = self.get('tag')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        structure_name = self.get('structure_name', 'minimum_relaxed')
        out_name = self.get('out_name', 'PPES@minimum')
        
        d_list = self.get('distance_list')
        
        interface_dict = GetInterfaceFromDB(name, db_file, functional)
        
        calc_PPES = True
        if interface_dict.get('PPES') is not None:
            if interface_dict['PPES'].get(out_name) is not None:
                print('\n A PPES-object with out_name: '+out_name+
                      '\n has already been created in the interface entry: '+
                      name+'\n for the '+functional+' functional.')
                calc_PPES = False
                
        if calc_PPES:
            SWF = CalcPPES_SWF(interface_name = name,
                               functional = functional,
                               distance_list = d_list,
                               out_name = out_name,
                               structure_name = structure_name,
                               spec = fw_spec)
            return FWAction(additions=SWF, update_spec=fw_spec)
        else:
            return FWAction(update_spec=fw_spec)

@explicit_serialize
class FT_DoPPESCalcs(FiretaskBase):
    """Start static PPES calculations in parallel for each given distance.
    
    First part of the PPES subworkflow (followed by FT_FitPPES). Takes the
    interface from the high-level database, moves that slabs to the distances
    passed in the distance_list and submits the calculations in parallel.
    
    Parameters
    ----------
    interface_name : str
        Name of the interface in the high-level database.
    functional : str
        Which functional to use; has to be 'PBE' or 'SCAN'.
    tag : str
        Unique tag to identify the calculations.
    distance_list : list of float, optional
        Modification of the equilibrium distance between the slabs.
    db_file : str
        Full path to the db.json file that should be used. Defaults to
        '>>db_file<<', to use env_chk.
    structure_name : str, optional
        Name of the structure in the interface entry to the high-level database
        for which the PPES should be calculated. The default is
        'minimum_relaxed'.
    out_name : str, optional
        Name for the PPES data in the high-level database. The default is
        'PPES@minimum'.
        
    Returns
    -------
    FWActions that produce a detour workflow with static PPES calculations.
    """
    required_params = ['interface_name', 'functional', 'tag', 'distance_list']
    optional_params = ['db_file', 'structure_name', 'out_name']
    def run_task(self, fw_spec):
        name = self.get('interface_name')
        functional = self.get('functional')
        tag = self.get('tag')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        structure_name = self.get('structure_name', 'relaxed_structure@min')
        
        d_list = self.get('distance_list')
        
        interface_dict = GetInterfaceFromDB(name, db_file, functional)            
        
        comp_params = interface_dict['comp_parameters']
        struct = Structure.from_dict(interface_dict[structure_name])
        sites_to_shift = []
        for i, s in enumerate(struct.sites):
            if s.c > 0:
                sites_to_shift.append(i)
        
        FW_list=[]
        for d in d_list:
            label = tag + '_PPES_' + str(d)
            #Make sure that there are no NoneTypes in the site_properties!
            struct_d = CleanUpSitePorperties(struct.copy())
            struct_d.translate_sites(indices=sites_to_shift,
                                     vector=[0,0,d],
                                     frac_coords=False, to_unit_cell=False)
            
            vis = GetCustomVaspStaticSettings(struct_d, comp_params,
                                              'slab_from_scratch')
                           
            FW = StaticFW(structure=struct_d, vasp_input_set=vis,
                                  name=label)
            FW_list.append(FW)
            
        
        WF = Workflow(FW_list, name='PPES calcs for: '+name)
        PPES_Calcs_WF = add_modify_incar(WF)
        
        return FWAction(detours = PPES_Calcs_WF)
    
@explicit_serialize 
class FT_FitPPES(FiretaskBase):
    """Fit PPES results with a UBER function and save them in the database.
    
    Second part of the PPES subworkflow (preceded by FT_DoPPESCalcs). Use
    a tag to read out the PPES calculation results from the database and then
    fit the results to an UBER function. Save everything in the high-level
    database.
    
    Parameters
    ----------
    interface_name : str
        Name of the interface in the high-level database.
    functional : str
        Which functional to use; has to be 'PBE' or 'SCAN'.
    tag : str
        Unique tag to identify the calculations.
    distance_list : list of float, optional
        Modification of the equilibrium distance between the slabs.
    db_file : str
        Full path to the db.json file that should be used. Defaults to
        '>>db_file<<', to use env_chk.
    out_name : str, optional
        Name for the PPES data in the high-level database. The default is
        'PPES@minimum'.
        
    """
    required_params = ['interface_name', 'functional', 'tag', 'distance_list']
    optional_params = ['db_file', 'out_name']
    def run_task(self, fw_spec):
        name = self.get('interface_name')
        functional = self.get('functional')
        tag = self.get('tag')
        db_file = self.get('db_file')
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        out_name = self.get('out_name', 'PPES@minimum')
        d_list = self.get('distance_list')
        
        
        DB = GetDB(db_file)
        
        d_E_array=[]
        for d in d_list:
            calc_label = tag + '_PPES_' + str(d)
            #Get energy from vasp run for this distance
            vasp_calc = DB.tasks.find_one({'task_label': calc_label})
            energy = vasp_calc['output']['energy']
            d_E_array.append([d, energy])
        
        tribo_db = GetHighLevelDB(db_file)
        coll = tribo_db[functional+'.interface_data']
        
        
        popt, perr = PPES_UBER(distance_energy_array=d_E_array)
        coll.update_one({'name': name},
                        {'$set': {'PPES':
                                    {out_name: 
                                    {'distance_Energy_array': d_E_array,
                                     'UBER': {'G': popt[0],
                                              'l': popt[1],
                                              'sigma_G': perr[0],
                                              'sigma_l': perr[1]}}}}})
        
                
        
            
def UBER(x, G, l):
    """
    Define the UBER function.
    
    Universal binding energy relation.

    Parameters
    ----------
    x : float or array of float
        distance(s).
    G : float
        binding energy at UBER minimum.
    l : float
        critcal length (location of inflection point).

    Returns
    -------
    float
        UBER binding energy at position(s) x.

    """
    return G * (1 -  (1 + x/l) * np.exp(-x/l))

def PPES_UBER(distance_energy_array):
    """
    Fit PPES data to an UBER relation.
    
    Takes a energy vs distance array and fits and UBER function to the data.

    Parameters
    ----------
    distance_energy_array : 2-D list of floats.
        Distance and energy pairs.

    Returns
    -------
    popt : list of floats
        Results of the fit to the UBER function (popt[0]=G, popt[1]=l).
    perr : float
        One standard deviation error.

    """
    #Fit UBER only to attractive part
    d_E_array=np.asarray(distance_energy_array)
    Fit_data = d_E_array[d_E_array[:,0] >= 0.0,:]
    z = Fit_data[:,0]
    E = Fit_data[:,1]-min(Fit_data[:,1])
    Initial_guess = [-min(E),0.75]
    popt, pcov = curve_fit(UBER, z, E, p0=Initial_guess)
    perr = np.sqrt(np.diag(pcov))
    return popt, perr