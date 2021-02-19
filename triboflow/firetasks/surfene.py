#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:54:53 2021

Firetasks to calculate the surface energy for a bulk and a given orientation.

@author: glosi000
"""

__author__ = 'Gabriele Losi'
__copyright__ = 'Prof. M.C. Righi, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 2nd, 2021'

import os

import numpy as np
from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks import FiretaskBase, FWAction

from triboflow.utils.database import Navigator
from triboflow.firetasks.slabs import read_runtask_params, get_multiple_info_from_struct_dict, write_multiple_dict_for_db
from triboflow.utils.errors import SurfaceEnergyError

currentdir = os.path.dirname(__file__)

# ============================================================================
# Firetasks
# ============================================================================

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


# ============================================================================
# Functions
# ============================================================================

def calculate_surface_energy(output_list, sym_surface=True):
    """
    Calculate the surface energy by passing a list containing all the dictionary
    with the output data. The first element is treated to be the bulk oriented
    along a specific miller index direction.

    Parameters
    ----------
    output_list : list of dict
        List of dictionaries where each one contains the output of a VASP
        simulation. Basically output_list can be whatever dictionary, fundamental
        keys that should be present are however: structure, energy, energy_per_atom.
    
    sym_surface : bool, optional
        If the surfaces are to be considered symmetric or not.
    
    Returns
    -------
    surfene : list of floats
        List containing the surface energies calculated from output_list.
    """

    # Take out the energy per atom of the bulk
    energy_per_atom = output_list[0]['energy_per_atom']

    # Loop over the slab elements of output_list and calculate surface energy
    surfene = np.array([])
    for d in output_list[1:]:
        energy = d['energy']
        nsites = d['nsites']
        surfene.append(energy - energy_per_atom * nsites)

    # Divide the surface energies by two if the surfaces are symmetric
    if sym_surface:
        surfene /= 2.

    return surfene
