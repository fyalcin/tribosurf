#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""

from fireworks import LaunchPad, Firework, Workflow
from fireworks.core.rocket_launcher import rapidfire
from triboflow.workflows.main import Heterogeneous_WF


inputs = {'material_1': {'formula': 'C',
                         'miller': '001',
                         'mp_id': 'mp-66',
                         'min_vacuum': 25,
                         'min_thickness': 10
                         },
          'material_2': {'formula': 'Au',
                         'miller': '111',
                         'min_vacuum': 25,
                         'min_thickness': 10
                         },
          'computational_params':{'functional': 'SCAN',
                                  'energy_tolerance': 0.001,
                                  'volume_tolerance': 0.001,
                                  'BM_tolerance': 0.01,
                                  'use_vdw': 'False'},
          'interface_params':{'interface_distance': 2.5,
                              'max_area': 500,
                              'r1r2_tol': 0.05
                              }
          }



WF = Heterogeneous_WF(inputs)

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
#rapidfire(lpad)
