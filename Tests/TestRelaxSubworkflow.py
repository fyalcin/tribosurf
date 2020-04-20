#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Common Workflows to be used by the FireFlow Project
Created on Wed Apr 15 17:04:35 2020

@author: mwo
"""

from fireworks import LaunchPad
from fireworks.core.rocket_launcher import rapidfire
from CommonWorkflows import Relax_SWF
from HelperFunctions import GetLowEnergyStructure, GetGapFromMP

struct , mp_id = GetLowEnergyStructure('NiO')

prev_spec = {'some_data': {'a': 1, 'b':2},
             'KPOINTS': 'previous_KPTS',
             'relax_loc': './RELAX_THIS_THING'}

bandgap = GetGapFromMP(mp_id)
if bandgap > 0.2:
    is_metal = False
else:
    is_metal = True

parameters = {'use_vdw': True,
              'functional': 'SCAN',
              'k_distance': 35,
              'encut': 350,
              'is_metal': is_metal}

wf = Relax_SWF(struct, parameters, 'bulk_pos_relax', prev_spec)
#fw = OptimizeFW(struct, override_default_vasp_params=uis)
#fw2 = RelaxFW(name='print spec')

lpad = LaunchPad.auto_load() # loads this based on the FireWorks configuration
#lpad.add_wf(Workflow([fw, fw2], {fw: [fw2]}))
lpad.add_wf(wf)


rapidfire(lpad)