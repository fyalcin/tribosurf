#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:54:53 2021

Firetasks to calculate the surface energy for a slab along a given orientation.

The module contains the following Firetasks:

** Surface Energy evaluation **
    - FT_SurfaceEnergy
    Calculate the surface energy between of one or more slabs with respect to
    the bulk structure having the same orientation.

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna

"""

__author__ = 'Gabriele Losi'
__copyright__ = 'Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 2nd, 2021'

import os

import numpy as np
from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks import FiretaskBase, FWAction

from triboflow.phys.surface_energy import calculate_surface_energy

from triboflow.utils.database import Navigator
from triboflow.firetasks.slabs import read_runtask_params, get_multiple_info_from_struct_dict, write_multiple_dict_for_db
from triboflow.utils.errors import SurfaceEnergyError

currentdir = os.path.dirname(__file__)


@explicit_serialize
class FT_SurfaceEnergy(FiretaskBase):
    """
    Calculate the surface energy between a bulk structure and one or more of the
    slabs oriented in the same way.
    
    """
    
    _fw_name = 'Surface Energy calculation'
    required_params = ['mp_id', 'collection', 'miller', 'entry']
    optional_params = ['db_file', 'database']
    
    def run_task(self, fw_spec):
        """ Run the Firetask.
        """

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_runtask_params(self, fw_spec, self.required_params,
                                self.optional_params, default_file=dfl,
                                default_key="SurfaceEnergy")

        # Calculate the surface energies from the provided entries
        surfene = self.get_surfene(p)

        # Store surface energies in DB
        self.store_to_db(surfene, p)

        return FWAction(update_spec=fw_spec)        

        def get_surfene(self, p):
        
            self._check_errors(p)
                
            # Call the navigator for retrieving the dictionary out of the DB
            nav = Navigator(db_file=p['db_file'], high_level=p['database'])
            dic = nav.find_data(collection=p['collection'], 
                                filter={'mpid': p['mp_id'], 'miller': p['miller']})
            
            # Extract the output dictionary containing energies and get surfene
            output_list = get_multiple_info_from_struct_dict(dic, p['entry'])
            surfene = calculate_surface_energy(output_list, sym_surface=True)

            return surfene
        
        def store_to_db(self, surfene, p):

            # Extract the surface energies from the provided dictionary entry
            entry = []
            for n in p['entry'][1:]:
                entry.append(n.append('surface_energy'))
            
            # Prepare the list of dictionaries to be stored in the database
            info_dict = write_multiple_dict_for_db(surfene, entry)
        
            # Prepare the database and options where to store data
            nav = Navigator(db_file=p['db_file'], high_level=p['database_to'])
 
            # Effectively store the data
            for d in info_dict:
                nav.update_data(p['collection'], 
                                {'mpid': p['mp_id'], 'miller': p['miller']}, 
                                {'$set': d})

        def _check_errors(self, p):
            SurfaceEnergyError.check_collection(p['collection'])
