#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 14:54:53 2021

Firetasks to calculate the surface energy for a slab along a given orientation.

The module contains the following Firetasks:

    - FT_SurfaceEnergy
    Calculate the surface energy between of one or more slabs with respect to
    the bulk structure having the same orientation.

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna

"""

__author__ = 'Gabriele Losi'
__copyright__ = 'Copyright 2021, Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 2nd, 2021'

import os

from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks import FiretaskBase, FWAction

from triboflow.phys.surface_energy import calculate_surface_energy

from triboflow.utils.database import Navigator
from triboflow.utils.utils import (
    read_runtask_params,
    get_multiple_info_from_dict,
    write_multiple_dict
)
from triboflow.utils.errors import SurfaceEnergyError


currentdir = os.path.dirname(__file__)


@explicit_serialize
class FT_SurfaceEnergy(FiretaskBase):
    """
    Calculate the surface energy between an oriented bulk structure and one or
    more of its slabs, oriented at the same way. This Firetasks simply combines
    the energies of the structures mentioned above, which need to be already
    calculated. The DB field is identified by means of `db_file`, `database` and
    `collection`. The entry should locate at the nested level of the subdicts
    present in field, where there is 'energy_per_atom' (bulk), 'energy' (slab). 
    Argument `entry` must be a list of lists: the first one identifies
    the oriented bulk, the other ones point to slabs. Surface energies are then 
    calculated and stored at the same level of the entries. 
    Warning: All of these structure should be located in the same field.

    mp_id : str
        MP-ID of the structure from the MP database, used to identify the
        structures in local databases.

    collection : str
        Collection where the structure is present. The default collections 
        used in the "tribchem" database are identified as functional+string, 
        where string is usually: '.bulk_data', 'slab_data', 'interface_name'.

    miller : list of int
        Miller indexes (h, k, l) used to identify unequivocally a slab within
        the database.

    entry : list of lists
        Location of the energy data to be analyzed from database+collection.

    db_file : str or None, optional
        Path to the location of the database. If nothing is provided it will be
        searched by env_check from Atomate. The default is None.

    database : str, optional
        Name of the database where the structure will be retrieved. 
        The default is "tribchem".
    
    Examples
    --------
    You might have a field like that:

    {
        "_id" : ...,
        "mp_id" : "mp-126",
        "miller" : [1, 1, 1],
        ...,

        "oriented_bulk" : {
            "energy_per_atom" : 1
        }
        "slab_1" : {
            "energy" : 10
        }
        "slab_2" : {
            "data" : {
                "energy" : 8
            }
        }
    }

    With this situation, your input arguments should be:
        - mp_id : "mp-126"
        - miller : [1, 1, 1]
        - entry = [['oriented_bulk'], ['slab_1'], ['slab_2', 'data']]
    
    At the end of the job the field will be updated as:

        {
        "_id" : ...,
        "mp_id" : "mp-126",
        "miller" : [1, 1, 1],
        ...,

        "oriented_bulk" : {
            "energy_per_atom" : 1
        }
        "slab_1" : {
            "energy" : 10
        }
        "slab_2" : {
            "data" : {
                "energy" : 8
                "surface_energy"
            }
        }
    }

    """
    
    _fw_name = 'Surface Energy calculation'
    required_params = ['mp_id', 'collection', 'miller', 'entry']
    optional_params = ['db_file', 'database']
    
    def run_task(self, fw_spec):
        """ Run the Firetask.
        """

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/../defaults.json'
        p = read_runtask_params(self, fw_spec, self.required_params,
                                self.optional_params, default_file=dfl,
                                default_key="SurfaceEnergy")

        # Calculate the surface energies from the provided entries
        surfene = self.get_surfene(p)

        # Store surface energies in DB
        self.store_to_db(surfene, p)

        return FWAction(update_spec=fw_spec)        

    def get_surfene(self, p):
        """
        Extract the entry from the field, find the energies and calculate surfene.

        Parameters
        ----------
        p : dict
            All the input parameters of the Firetasks, placed in a dictionary.

        Returns
        -------
        surfene : list of floats
            List of surface energies calculated of the slabs.

        """
        
        # Check for errors in the length of entry
        SurfaceEnergyError.check_entry(p['entry'])
            
        # Call the navigator for retrieving the dictionary out of the DB
        nav = Navigator(db_file=p['db_file'], high_level=p['high_level_db'])
        dic = nav.find_data(collection=p['collection'], 
                            fltr={'mpid': p['mp_id'], 'miller': p['miller']})

        # Extract the output dictionary containing energies and get surfene
        try:
            output_list = get_multiple_info_from_dict(dic, p['entry'])
        except:
            raise SurfaceEnergyError('Some problems occurred in output list, '
                                     'probably results are not stored correctly '
                                     'in database: {}, collection: {}'
                                     .format(p['database'], p['collection']))
        SurfaceEnergyError.check_output(output_list)
        
        surfene = calculate_surface_energy(output_list, sym_surface=True)

        return surfene

    def store_to_db(self, surfene, p):
        """
        Store surface energies to DB.

        Parameters
        ----------
        surfene : list
            List of surface energies calculated of the slabs.

        p : dict
            All the input parameters of the Firetasks.
        """

        # Extract the surface energies from the provided dictionary entry
        entry = p['entry'][1:]
        _ = [n.append('surface_energy') for n in entry]  # Add nested key
        
        # Prepare the list of dictionaries to be stored in the database      
        info_dict = write_multiple_dict(surfene, entry)
    
        # Prepare the database and store the data one by one
        nav = Navigator(db_file=p['db_file'], high_level=p['high_level_db'])
        for d in info_dict:
            nav.update_data(p['collection'], 
                            {'mpid': p['mp_id'], 'miller': p['miller']}, 
                            {'$set': d})        
