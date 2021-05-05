#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""

from fireworks import LaunchPad, Firework, Workflow
from fireworks.core.rocket_launcher import rapidfire
from triboflow.workflows.main import heterogeneous_wf


inputs = {'material_1': {'formula': 'Al',
                         'miller': '111',
                         'mp_id': 'mp-134',
                         'vacuum': 25,
                         'thick_min': 3,
                         'thick_max': 12,
                         'thick_incr': 1
                         },
          'material_2': {'formula': 'Cu',
                         'miller': '111',
                         'mp_id': 'mp-30',
                         'vacuum': 25,
                         'thick_min': 3,
                         'thick_max': 12,
                         'thick_incr': 1
                         },
          'computational_params':{'functional': 'PBE',
                                  'volume_tolerance': 0.001,
                                  'BM_tolerance': 0.01,
                                  'use_vdw': 'False',
                                  'surfene_thr': 0.01},
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