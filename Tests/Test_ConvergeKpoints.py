#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
Created on Wed Apr 15 17:04:35 2020

@author: mwo
"""

from fireworks import LaunchPad, Firework, Workflow
from fireworks.core.rocket_launcher import rapidfire
from CommonWorkflows import ConvergeKpoints_SWF
from HelperFunctions import GetLowEnergyStructure, GetGapFromMP

# =============================================================================
# THIS IS A SIMPLE TEST AS A COMPLETE WORKFLOW:

struct , mp_id = GetLowEnergyStructure('Al')

prev_spec = {'some_data': {'a': 1, 'b':2},
              'KPOINTS': 'previous_KPTS',
              'relax_loc': './RELAX_THIS_THING'}

bandgap = GetGapFromMP(mp_id)
if bandgap > 0.1:
    is_metal = False
    k_distance_start = 10
else:
    is_metal = True
    k_distance_start = 40

parameters = {'use_vdw': False,
              'functional': 'PBE',
              #'k_distance': 25,
              'encut': 350,
              'is_metal': is_metal}

kinfo_loc = ['material_1','Kpoint_Convergence']
wf = ConvergeKpoints_SWF(struct, comp_parameters=parameters,
                         out_loc=['material_1'], k_dist_start=k_distance_start,
                         to_pass=['material_1'], k_dist_incr=5, spec=prev_spec)
#fw = OptimizeFW(struct, override_default_vasp_params=uis)
#fw2 = RelaxFW(name='print spec')

lpad = LaunchPad.auto_load() # loads this based on the FireWorks configuration
#lpad.add_wf(Workflow([fw, fw2], {fw: [fw2]}))
lpad.add_wf(wf)


rapidfire(lpad)
