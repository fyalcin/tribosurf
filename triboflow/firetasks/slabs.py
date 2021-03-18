#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 15:09:03 2021

Firetasks to generate and manipulate crystalline slabs.

The module contains the following Firetasks:

** Slab generation **

    - FT_GenerateSlabs
    Generate a list of slabs out of a structure provided as input and store
    them in a given location inside the database.

    Author: Gabriele Losi (glosi000)
    Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna

"""

__author__ = 'Gabriele Losi'
__copyright__ = 'Copyright 2021, Prof. M.C. Righi, TribChem, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 2nd, 2021'


import os

from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks import FiretaskBase, FWAction

from triboflow.phys.solid_state import generate_slabs
from triboflow.utils.database import Navigator
from triboflow.utils.utils import (
    read_runtask_params,
    write_one_dict
)
from triboflow.utils.errors import GenerateSlabsError


currentdir = os.path.dirname(__file__)


@explicit_serialize
class FT_GenerateSlabs(FiretaskBase):
    """
     Generate a slab or a list of slabs out of a given structure. Parameters 
     hat are taken into account to generate the possible different slabs are: 
     miller, thickness, vacuum, entry. The slabs are generated with 
     SlabGenerator and stored in the database.
    
    """

    required_params = ['structure', 'mp_id', 'miller', 'collection']
    optional_params = ['db_file', 'database', 'thickness', 'thick_bulk', 'vacuum',
                       'symmetrize', 'ext_index', 'in_unit_planes', 'entry']

    def run_task(self, fw_spec):
        """ Run the Firetask.
        """ 

        # Define the json file containing default values and read parameters
        dfl = currentdir + '/defaults_fw.json'
        p = read_runtask_params(self, fw_spec, self.required_params, self.optional_params,
                                default_file=dfl, default_key="GenerateSlabs")

        # Generate the slabs for each structure passed as input
        miller, thickness, vacuum, entry = self._convert_to_list(p)

        # Generate the slabs, taking into account miller, thickness, vacuum
        slabs = generate_slabs(structure=p['structure'], 
                               miller=miller, 
                               thickness=thickness, 
                               vacuum=vacuum, 
                               thick_bulk=p['thick_bulk'],
                               ext_index=p['ext_index'],
                               in_unit_planes=p['in_unit_planes'])
        
        # Store the slab structure in collection.field.entry in db, within db_file
        self.structure_in_db(slabs, miller, entry, p)
        
        return FWAction(update_spec=fw_spec)
    
    def _convert_to_list(self, p):

        GenerateSlabsError.check_inputs(p['miller'], p['thickness'], 
                                        p['vacuum'], p['entry'])
        
        miller = p['miller']
        if not all([isinstance(x, list) for x in miller]):
            miller = [miller]
        
        thickness = p['thickness']
        if not isinstance(thickness, list):
            thickness = [thickness] * len(miller)
        
        entry = p['entry']
        if not isinstance(entry, list):
            entry = [entry] * len(miller)

        vacuum = p['vacuum']
        if not isinstance(vacuum, list):
            vacuum = [vacuum] * len(miller)

        return miller, thickness, vacuum, entry
    
    def structure_in_db(self, slabs, miller, entry, p):

        nav = Navigator(p['db_file'], p['database'])
        
        # Store unrelaxed data in the Database
        for s, hkl, en in zip(slabs, miller, entry):
            # Clean the data and create a dictionary with the given path
            update_data = write_one_dict(s.as_dict(), en)
            nav.update_data(collection=p['collection'], 
                            fltr={'mpid': p['mp_id'], 'miller': hkl},
                            new_values={'$set': update_data},
                            upsert=True)
