#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 10:37:34 2020

@author: mwo
"""

from fireworks import LaunchPad
from fireworks.core.rocket_launcher import rapidfire
from CommonWorkflows import Heterogeneous_WF



#Create the input dictionary

inputs = {'material_1': {'formula': 'Al',
                         'miller': '111',
                         'min_vacuum': 25,
                         'min_thickness': 5
                         },
          'material_2': {'formula': 'Au',
                         #'mp_id': 'mp-13',
                         'miller': '110',
                         'min_vacuum': 25,
                         'min_thickness': 6
                         },
          'computational_params':{'functional': 'PBE',
                                  'energy_tolerance': 0.01,
                                  'volume_tolerance': 0.01,
                                  'BM_tolerance': 0.01,
                                  'use_vdw': 'False'},
          'interface_params':{'interface_distance': 1.5,
                              'max_area': 500,
                              'r1r2_tol': 0.05
                              }
          }


lpad = LaunchPad.auto_load() # loads this based on the FireWorks configuration
lpad.add_wf(Heterogeneous_WF(inputs))

#rapidfire(lpad)