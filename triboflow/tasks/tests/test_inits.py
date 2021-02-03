#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 09:14:37 2021

Test the initializing Firetasks for setting up a Workflow.

@author: glosi000
"""

from triboflow.tasks.init_check import FTCheckInput
from triboflow.tasks.init_fws import InitFWs
from fireworks import LaunchPad, Workflow
from fireworks.core.rocket_launcher import rapidfire

inputs = {'material_1': {'formula': 'Pt',
                         'miller': '111',
                         'min_vacuum': 35,
                         'min_thickness': 10
                         },
          'material_2': {'formula': 'Ag',
                         'miller': '111',
                         'mpid': 'mp-124',
                         'min_vacuum': 35,
                         'min_thickness': 10
                         },
          'computational_params':{'functional': 'PBE',
                                  'energy_tolerance': 0.001,
                                  'volume_tolerance': 0.001,
                                  'BM_tolerance': 0.01,
                                  'use_vdw': 'False'},
          'interface_params':{'interface_distance': 2.5,
                              'max_area': 500,
                              'r1r2_tol': 0.05
                              }
          }

mat_1 = inputs.get('material_1')
mat_2 = inputs.get('material_2')
comp_params = inputs.get('computational_params')
inter_params = inputs.get('interface_params')

# A little bit slow, recover data from online MP Database
#structs, mpids, functional = materials_from_mp(inputs)

#fw = InitFWs.checkinp_hetero_interface(mat_1, mat_2, comp_params, inter_params)
fw = InitFWs.checkinp_homo_interface(mat_1, comp_params, inter_params)
wf = Workflow([fw], name='test_init')

lpad = LaunchPad.auto_load()
lpad.add_wf(wf)
rapidfire(lpad)
