#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
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
if bandgap > 0.1:
    is_metal = False
else:
    is_metal = True

parameters = {'use_vdw': False,
              'functional': 'PBE',
              #'k_distance': 25,
              #'encut': 350,
              'is_metal': is_metal}

Relaxed_structure_loc = ['structures',
                         'relaxed_structures',
                         struct.composition.reduced_formula]
wf = Relax_SWF(struct, parameters, 'bulk_pos_relax', Relaxed_structure_loc,
               prev_spec)
#fw = OptimizeFW(struct, override_default_vasp_params=uis)
#fw2 = RelaxFW(name='print spec')

lpad = LaunchPad.auto_load() # loads this based on the FireWorks configuration
#lpad.add_wf(Workflow([fw, fw2], {fw: [fw2]}))
lpad.add_wf(wf)


#rapidfire(lpad)