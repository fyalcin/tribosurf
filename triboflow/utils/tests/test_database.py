#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 15:15:11 2021

Classes to manage data from local and online DataBases at a high level.

@author: omarchehaimi
"""

__author__ = 'Omar Chehaimi'
__copyright__ = 'Prof. M.C. Righi, TribChem, ERC-SLIDE, University of Bologna'
__contact__ = 'clelia.righi@unibo.it'
__date__ = 'February 1st, 2021'

import os
from pathlib import Path, PurePosixPath
from datetime import datetime

from triboflow.utils.database import Navigator, NavigatorMP

# PurePosixPath gets the first level parten directory
project_folder = os.getcwd()
db_folder_obj = PurePosixPath(project_folder)
db_file = str(db_folder_obj.parent) + '/config/db.json'

# Retrive first some data from Materials Project database and do some tests
nav_mp = NavigatorMP()
mp_id = nav_mp.get_mpid_from_formula('NaCl')
print('NaCl mp id: ', mp_id)
les = nav_mp.get_low_energy_structure('NaCl')
print('Low energy structure: ', les)
energy_properties = nav_mp.get_property_from_mp(
    mp_id,
    ['energy', 'energy_per_atom'])
print('Energy: ', energy_properties['energy'], 
      'Energy per atom: ', energy_properties['energy_per_atom'])

# Test the writing/reading operations in the database
nav = Navigator()
# Insert
nav.insert_data('Collection_test', 
                {'tag': 0,
                 'chem_formula': 'NaCl',
                 'Energy': 10,
                 'Energy per atom': 10
                })
# Insert again the same data and check if the method find_data works.
# It should print that the data already exist and not save the data.
nav.insert_data('Collection_test', 
                {'tag': 0,
                 'chem_formula': 'NaCl',
                 'Energy': 10,
                 'Energy per atom': 10
                })
# Insert again the same data and check if the method find_data works.
nav.insert_data('Collection_test', 
                {'tag': 0,
                 'chem_formula': 'NaCl',
                 'Energy': 10,
                 'Energy per atom': 10
                }, duplicates=True)
# Insert the same data twice for testing if does not save duplicates data
nav.insert_many_data('Collection_test', 
                     [{'tag': 0,
                       'chem_formula': 'NaCl',
                       'Energy': 10,
                       'Energy per atom': 10
                      }, {'tag': 0,
                       'chem_formula': 'NaCl',
                       'Energy': 10,
                       'Energy per atom': 10
                      }])
# Insert the same data twice for testing delete many data later
nav.insert_many_data('Collection_test', 
                     [{'tag': 0,
                       'chem_formula': 'NaCl',
                       'Energy': 10,
                       'Energy per atom': 10
                      }, {'tag': 0,
                       'chem_formula': 'NaCl',
                       'Energy': 10,
                       'Energy per atom': 10
                      }], duplicates=True)
# Update one
nav.update_data('Collection_test',
                {'tag': 0,
                 'chem_formula': 'NaCl',
                 'Energy': 10,
                 'Energy per atom': 10
                },
                {'$set': {'tag': 0,
                 'chem_formula': 'NaCl',
                 'Energy': 1000,
                 'Energy per atom': 10
                }})
# Update many
nav.update_many_data('Collection_test',
                     {'tag': 0,
                      'chem_formula': 'NaCl',
                      'Energy': 1000,
                      'Energy per atom': 10
                      },
                     {'$set': {'tag': 0,
                      'chem_formula': 'NaCl',
                      'Energy': 10,
                      'Energy per atom': 10
                     }})   
# Test find many data
data = nav.find_many_data('Collection_test', 
                          {'tag': 0,
                          'chem_formula': 'NaCl',
                          'Energy': 10,
                          'Energy per atom': 10
                          })
print('Retrived data from the database: ', data)
# Test find many data
data = nav.find_data('Collection_test', 
                     {'tag': 0,
                      'chem_formula': 'NaCl',
                      'Energy': 10,
                      'Energy per atom': 10
                      })
print('Retrived data from the database: ', data)
# Remove the data
nav.delete_many_data('Collection_test', 
                     {'tag': 0,
                      'chem_formula': 'NaCl',
                      'Energy': 10,
                      'Energy per atom': 10
                     })
# Insert data for removing with delete one
nav.insert_data('Collection_test', 
                {'tag': 0,
                 'chem_formula': 'NaCl',
                 'Energy': 20,
                 'Energy per atom': 10
                })
# Removing the previous added data
nav.delete_data('Collection_test', 
                {'tag': 0,
                 'chem_formula': 'NaCl',
                 'Energy': 20,
                 'Energy per atom': 10
                })
# Try to the data which have been just removed
nav.find_data('Collection_test',
              {'tag': 0,
               'chem_formula': 'NaCl',
               'Energy': 20,
               'Energy per atom': 10
              })
# Insert the same data twice for testing drop data later
nav.find_many_data('Collection_test', 
                   {'tag': 0,
                    'chem_formula': 'NaCl',
                    'Energy': 10,
                    'Energy per atom': 10
                   })
# Insert the same data twice for testing drop data later
nav.insert_many_data('Collection_test', 
                     [{'tag': 0,
                       'chem_formula': 'NaCl',
                       'Energy': 10,
                       'Energy per atom': 10
                      }, {'tag': 0,
                          'chem_formula': 'NaCl',
                          'Energy': 10,
                          'Energy per atom': 10
                      }])
# Drop all data
nav.drop_data('Collection_test')
