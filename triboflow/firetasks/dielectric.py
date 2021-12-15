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
class FT_StartDielectric(FiretaskBase):
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
        from triboflow.workflows.subworkflows import dielectric_constant_swf
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