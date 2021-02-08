#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""

from fireworks import LaunchPad, Firework, Workflow
from fireworks.core.rocket_launcher import rapidfire
from triboflow.workflows.main import heterogeneous_wf


inputs = {'material_1': {'formula': 'MgO',
                         'miller': '001',
                         'mp_id': 'mp-1918',
                         'min_vacuum': 25,
                         'min_thickness': 6
                         },
          'material_2': {'formula': 'FeRh',
                         'miller': '001',
                         'mp_id': 'mp-1265',
                         'min_vacuum': 25,
                         'min_thickness': 6
                         },
          'computational_params':{'functional': 'PBE',
                                  'energy_tolerance': 0.001,
                                  'volume_tolerance': 0.001,
                                  'BM_tolerance': 0.01,
                                  'use_vdw': 'False'},
          'interface_params':{'interface_distance': 2.5,
                              'max_area': 50,
                              'r1r2_tol': 0.2,
                              'max_mismatch': 0.1
                              }
          }



WF = heterogeneous_wf(inputs)

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
#rapidfire(lpad)