#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""

from fireworks import LaunchPad
from fireworks.core.rocket_launcher import rapidfire
from triboflow.workflows.main import heterogeneous_wf


inputs = {'material_1': {'formula': 'Al',
                         'miller': '111',
                         'mpid': 'mp-134',
                         'thick_min': 3,
                         'thick_max': 12,
                         'thick_incr': 1,
                         },
          'material_2': {'formula': 'C',
                         'miller': '111',
                         'mpid': 'mp-66',
                         'thick_min': 6,
                         'thick_max': 18,
                         'thick_incr': 2
                         },
          'computational_params':{'functional': 'PBE',
                                  'volume_tolerance': 0.001,
                                  'BM_tolerance': 0.01,
                                  'use_vdw': 'No',
                                  'surfene_thr': 0.01,
                                  'vacuum': 15},
          'interface_params':{'interface_distance': 'auto',
                              'max_area': 200.0,
                              'max_mismatch': 0.02
                              }
          }

WF = heterogeneous_wf(inputs)

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
#rapidfire(lpad)
