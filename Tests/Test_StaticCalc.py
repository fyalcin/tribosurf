#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 16:07:23 2020

@author: mwo
"""


from fireworks import LaunchPad, Firework, Workflow
from fireworks.core.rocket_launcher import rapidfire
from atomate.vasp.fireworks.core import StaticFW
from CommonFiretasks import FT_PrintSpec, FT_PassSpec, FT_GetEnergyFromDB
from CommonFireworks import ParsePrevVaspCalc_FW
from HelperFunctions import GetLowEnergyStructure


struct, mp_id = GetLowEnergyStructure('Al')

# =============================================================================
# FW_s = StaticFW(structure=struct)
# FW1 = Firework([FT_GetEnergyFromDB(label='static', formula='Al', out_loc=['results',
#                                                                    'Al',
#                                                                    'bulk',
#                                                                    'energy']),
#                FT_PrintSpec()])
# 
# WF = Workflow([FW_s, FW1], {FW_s: [FW1]})
# =============================================================================

calc_name = 'Al_bulk_static'
FW_s = StaticFW(structure=struct, name=calc_name)
FW_p = ParsePrevVaspCalc_FW(calc_name, ['results', calc_name], ['_all'])
FW_end = Firework(FT_PrintSpec())

WF = Workflow([FW_s, FW_p, FW_end], {FW_s: [FW_p], FW_p: [FW_end]},
              name = 'TestWF for static Calc')

lpad = LaunchPad.auto_load() # loads this based on the FireWorks configuration
lpad.add_wf(WF)

rapidfire(lpad)