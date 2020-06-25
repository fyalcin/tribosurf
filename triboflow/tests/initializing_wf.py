#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""

from fireworks import LaunchPad, Firework, Workflow
from fireworks.core.rocket_launcher import rapidfire
from triboflow.fireworks.common import CheckInputsFW


inputs = {'material_1': {'formula': 'C',
                         'miller': '001',
                         'mp_id': 'mp-66',
                         'min_vacuum': 25,
                         'min_thickness': 5
                         },
          'material_2': {'formula': 'Al',
                         #'mp_id': 'mp-13',
                         'miller': '111',
                         'min_vacuum': 25,
                         'min_thickness': 6
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

Start_WF_FW = CheckInputsFW(mat1_params = inputs['material_1'],
                            mat2_params = inputs['material_2'],
                            compparams = inputs['computational_params'],
                            interface_params = inputs['interface_params'],
                            FW_name = 'Check input parameters FW')

WF = Workflow.from_Firework(Start_WF_FW, name='Test Workflow')

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
rapidfire(lpad)
