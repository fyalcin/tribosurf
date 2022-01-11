#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""

from fireworks import LaunchPad
from fireworks.core.rocket_launcher import rapidfire
from triboflow.workflows.main import heterogeneous_wf


inputs = {'material_1': {'formula': 'ZnCu',
                         'miller': '111',
                         'mpid': 'mp-987',
                         'thick_min': 4,
                         'thick_max': 12,
                         'thick_incr': 1,
                         },
          'material_2': {'formula': 'WC',
                         'miller': [1,1,-2],
                         'mpid': 'mp-1894',
                         'thick_min': 4,
                         'thick_max': 12,
                         'thick_incr': 1
                         },
          'computational_params':{'functional': 'PBE',
                                  'volume_tolerance': 0.001,
                                  'BM_tolerance': 0.01,
                                  'use_vdw': 'No',
                                  'surfene_thr': 0.01,
                                  'vacuum': 12},
          'interface_params':{'max_area': 100,
                              'r1r2_tol': 0.1,
                              'max_mismatch': 0.05
                              }
          }

WF = heterogeneous_wf(inputs)

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
#rapidfire(lpad)
