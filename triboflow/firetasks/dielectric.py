#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 16:11:53 2021

@author: mwo
"""

from fireworks import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from triboflow.utils.database import Navigator


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
        bandgap = output_data['output']['bandgap']
        
        # If we are dealing with a metal, set epsilon very high since it should be inf.
        if bandgap < 0.05:
            eps_average = 1000000
        
        return FWAction(update_spec={'epsilon_tensor': eps_tensor,
                                     'epsilon': eps_average})