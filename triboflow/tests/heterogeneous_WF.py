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
          'material_2': {'formula': 'Cu',
                         'miller': '111',
                         'mpid': 'mp-30',
                         'thick_min': 3,
                         'thick_max': 12,
                         'thick_incr': 1
                         },
          'computational_params':{'functional': 'PBE',
                                  'volume_tolerance': 0.01,
                                  'BM_tolerance': 0.01,
                                  'use_vdw': 'dftd2',
                                  'use_spin': 'Yes',
                                  'surfene_thr': 0.01,
                                  'vacuum': 20},
          'interface_params':{'max_area': 100,
                              'max_area_ratio_tol': 0.1,
                              'max_length_tol': 0.05
                              }
          }

WF = heterogeneous_wf(inputs)

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
#rapidfire(lpad)
