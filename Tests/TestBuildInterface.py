#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 16:02:54 2020

@author: mwo
"""

from fireworks import LaunchPad
from CommonFireworks import FW_MakeHeteroInterfaceFromScratch
from fireworks.core.rocket_launcher import rapidfire

launchpad = LaunchPad.auto_load()

Material_1 = {'formula': 'CuO',
              'miller': [0, 0, 1],
              'thickness': 7.0}

Material_2 = {'formula': 'GaAs',
              'miller': [0, 0, 1],
              'thickness': 7.0,
              #'MP_ID': 'mp-1178232',
              'vacuum': 25.0}

parameters = {'interface_distance': 1.5,
              'max_area': 100,
              'max_mismatch': 0.025,
              'max_angle': 2.5,
              'r1r2_tol': 0.025}

firework = FW_MakeHeteroInterfaceFromScratch(Material_1, Material_2, parameters)

launchpad.add_wf(firework)

rapidfire(launchpad)
