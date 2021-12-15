#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 16:11:53 2021

@author: mwo
"""
from uuid import uuid4

from pymatgen.core.structure import Structure

from fireworks import FWAction, FiretaskBase, Firework, Workflow
from fireworks.utilities.fw_utilities import explicit_serialize

from atomate.vasp.fireworks.core import StaticFW
from atomate.utils.utils import env_chk

from triboflow.utils.database import StructureNavigator, Navigator
from triboflow.utils.vasp_tools import get_custom_vasp_static_settings


@explicit_serialize
class FT_StartConvo(FiretaskBase):
    """ Starts a dielectric subworkflow.

    Starts a subworkflow that calculates and updates the dielectric constant
    for a given bulk material. Will also by default update all slabs with the
    same mpid.

    Parameters
    ----------
    mp_id : str
        MaterialsProject ID number for the material
    functional : str
        Functional with which the workflow is run. PBE or SCAN.
    db_file : str, optional
        Full path of the db.json file to be used. The default is to use
        env_chk to find the file.
    high_level_db : str or True, optional
        Name of the high_level database to use. Defaults to 'True', in which
        case it is read from the db.json file.
    update_bulk : bool, optional
        If the bulk entry for the given mpid should be updated.
        The default is True.
    update_slabs : bool, optional
        If the slab entries matching a given mpid should be updated (all miller
        indices. The default is False.
    """

    _fw_name = 'Start Encut or Kdensity Convergence'
    required_params = ['mp_id', 'functional']
    optional_params = ['db_file', 'high_level_db', 'update_bulk', 'update_slabs']

    def run_task(self, fw_spec):

        mp_id = self.get('mp_id')
        functional = self.get('functional')
        db_file = self.get('db_file')
        update_bulk = self.get('update_bulk', True)
        update_slabs = self.get('update_slabs', True)
        if not db_file:
            db_file = env_chk('>>db_file<<', fw_spec)
        hl_db = self.get('high_level_db', True)
            

        nav_structure = StructureNavigator(
            db_file=db_file, 
            high_level=hl_db)
        data = nav_structure.get_bulk_from_db(
            mp_id=mp_id, 
            functional=functional)
        
        structure_dict = data.get('structure_equiVol')
        if not structure_dict:
            structure_dict = data.get('primitive_structure')
            if not structure_dict:
                structure_dict = data.get('structure_fromMP')
                if not structure_dict:
                    raise LookupError('No structure found that can be used '
                                      'as input for the convergence swf.')
        structure = Structure.from_dict(structure_dict)
        comp_params = data.get('comp_parameters', {})
        flag = 'test_dielectric_WF_'+str(uuid4())
        SWF = dielectric_constant_swf(structure=structure,
                                      mpid=mp_id,
                                      flag=flag, 
                                      comp_parameters=comp_params,
                                      functional='PBE',
                                      db_file=db_file,
                                      hl_db=hl_db,
                                      update_bulk=update_bulk,
                                      update_slabs=update_slabs)
        return FWAction(detours=SWF, update_spec=fw_spec)


@explicit_serialize
class FT_GetEpsilon(FiretaskBase):
    """
    Firetask to read out the dielectric tensor from a vasp run.
    
    Parameters
    ----------
    label : str
        Task label to query for in the tasks collection of the Fireworks database.
    db_file : str, optional
        Full path of the db.json. The default is 'auto'.
    
    Returns
    -------
    FWAction that updates the spec with the dielectric tensor and the average
    value of its diagonal elements.
    
    """
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
    """
    Firetask to update computational parameters for bulk and/or slabs
    
    Parameters
    ----------
    mpid : str
        Material Project's material identifier ID.
    functional : str
        Functional for the identification of the high_level db.
    new_params : list
        List of strings that identify the new keys that should be written to
        the computational parameters.
    db_file : str, optional
        Full path of the db.json. The default is 'auto'.
    update_bulk : bool, optional
        If the bulk entry for the given mpid should be updated.
        The default is True.
    update_slabs : bool, optional
        If the slab entries matching a given mpid should be updated (all miller
        indices. The default is False.
    high_level_db : str or bool, optional
        If a string is given, the high-level database will be chosen based on
        that string. If True, the db.json file will be used to determine the
        name of the high_level_db. The default is True.

    """
    _fw_name = 'Update computational parameters'
    required_params = ['mpid', 'functional', 'new_params']
    optional_params = ['db_file', 'update_bulk', 'update_slabs',
                       'high_level_db']
    def run_task(self, fw_spec):

        mpid = self.get('mpid')
        functional = self.get('functional')
        new_params = self.get('new_params')
        update_bulk = self.get('update_bulk', True)
        update_slabs = self.get('update_slabs', False)
        db_file = self.get('db_file', 'auto')
        hl_db = self.get('high_level_db', True)
        
        nav_high = StructureNavigator(db_file=db_file, 
                                           high_level=hl_db)
        #get values from spec:
        new_data = {'$set': {}}
        for param in new_params:
            new_data['$set'][f'comp_parameters.{param}'] = fw_spec.get(param, None)
            
        if update_bulk:
            nav_high.update_data(collection=functional+'.bulk_data',
                                 fltr={'mpid': mpid},
                                 new_values=new_data,
                                 upsert=False)
           
        if update_slabs:
            nav_high.update_many_data(collection=functional+'.slab_data',
                                      fltr={'mpid': mpid},
                                      new_values=new_data,
                                      upsert=False)
        return
        

def dielectric_constant_swf(structure,
                            mpid,
                            flag, 
                            comp_parameters={}, 
                            spec={},
                            functional='PBE', 
                            db_file='auto',
                            hl_db=True,
                            update_bulk=True,
                            update_slabs=False):
    """
    Subworkflow that calculates dielectric properties and updates the comp_parameters.

    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        The structure for which the dielectric porperties are calculated.
    mpid : str
        Material Project's material identifier ID of the structure passed.
    flag : str
        Will be used to name the FW of the vasp calc and later to find it in the
        tasks collection.
    comp_parameters : dict, optional
        Dictionary of computational parameters. Used to set up the vasp input
        set. The default is {}.
    spec : dict, optional
        fw_spec that can be passed to the FWs. The default is {}.
    functional : str, optional
        Functional for the calculation. Usually SCAN or PBE. Used to select
        the output collection in the high level db. The default is 'PBE'.
    db_file : str, optional
        Path to a db.json file. If 'auto', the standard config folder is used.
        The default is 'auto'.
    hl_db : str or bool, optional
        If a string is given, the high-level database will be chosen based on
        that string. If True, the db.json file will be used to determine the
        name of the high_level_db. The default is True.
    update_bulk : bool, optional
        If the bulk entry for the given mpid should be updated.
        The default is True.
    update_slabs : bool, optional
        If the slab entries matching a given mpid should be updated (all miller
        indices. The default is False.

    Returns
    -------
    WF : fireworks.core.firework.Workflow
        The dielectric subworkflow.

    """
    
    #check if epsilon is already calculated for that material.
    try:
        nav_high = StructureNavigator(db_file=db_file, 
                                      high_level=hl_db)
        bulk_data = nav_high.get_bulk_from_db(mpid, functional)
        epsilon = bulk_data['comp_parameters'].get('epsilon', False)
    except:
        epsilon = False
    
    formula = structure.composition.reduced_formula
    wf_name = f'Dielectric calculation WF for {formula} {mpid}'
    
    if not epsilon:
        vis = get_custom_vasp_static_settings(structure,
                                              comp_parameters,
                                              'bulk_epsilon_from_scratch')
        Calc_Eps_FW = StaticFW(structure=structure,
                               name=flag,
                               vasp_input_set=vis)
    
        Get_Eps_FT = FT_GetEpsilon(label=flag, db_file=db_file)

        Update_Data_FT = FT_UpdateCompParams(mpid=mpid,
                                             functional=functional,
                                             new_params=['epsilon'],
                                             update_bulk=update_bulk,
                                             update_slabs=update_slabs,
                                             db_file=db_file,
                                             high_level_db=hl_db)
        Update_FW = Firework(tasks=[Get_Eps_FT, Update_Data_FT],
                             spec=spec,
                             name=flag+'_update_high_level')
    
        WF = Workflow([Calc_Eps_FW, Update_FW], {Calc_Eps_FW: [Update_FW]},
                      name=wf_name)
    else:
        spec.update({'epsilon': epsilon})
        Update_Data_FT = FT_UpdateCompParams(mpid=mpid,
                                             functional=functional,
                                             new_params=['epsilon'],
                                             update_bulk=update_bulk,
                                             update_slabs=update_slabs,
                                             db_file=db_file,
                                             high_level_db=hl_db)
        Update_FW = Firework(tasks=[Update_Data_FT],
                             spec=spec,
                             name=flag+'_update_high_level')
        WF = Workflow([Update_FW], name=wf_name)
    return WF