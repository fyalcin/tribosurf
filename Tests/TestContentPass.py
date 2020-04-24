#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 16:07:23 2020

@author: mwo
"""


from fireworks import LaunchPad, Firework, Workflow
from fireworks.core.firework import FWAction, FiretaskBase, Firework
from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks.core.rocket_launcher import rapidfire
from CommonFireworks import TestFW
from CommonFiretasks import FT_PrintSpec, FT_PassSpec, FT_CopyInSpec, FT_DoNothing

inputs = {'a': {'a1': 1},
          'b': {'a1': 2},
          'c': {'a': 'starta'},
          'd': {'a': 'startb'}
          }

# FW1 = Firework(FT_PrintSpec(), inputs)
# FW2 = Firework(FT_PassSpec(key_list = ['c']))
# FW3 = Firework(FT_PassSpec(key_list = ['d']))
# FW4 = Firework(FT_PrintSpec())

# Dependencies={FW1: [FW2, FW3],
#               FW2: [FW4],
#               FW3: [FW4]
#               }

# WF = Workflow([FW1, FW2, FW3, FW4], Dependencies)

FW1 = Firework(FT_PrintSpec(), inputs)
# FW2 = TestFW(in_loc=['a', 'a1'], out_loc=['c', 'b'], pass_list=['c'])
# FW3 = TestFW(in_loc=['b', 'a1'], out_loc=['d', 'b'], pass_list=['d'])
FW2 = Firework([FT_CopyInSpec(in_loc=['a', 'a1'], out_loc=['c', 'new']), 
                FT_DoNothing(), FT_PassSpec(key_list=['c'])])
#FW2a = Firework(FT_PassSpec(key_list=['c', 'a']))
FW3 = Firework([FT_CopyInSpec(in_loc=['b', 'a1'], out_loc=['d', 'new']), 
                FT_DoNothing(), FT_PassSpec(key_list=['d'])])
#FW3a = Firework(FT_PassSpec(key_list=['d']))
FW4 = Firework(FT_PrintSpec())

Dependencies={FW1: [FW2, FW3],
              FW2: [FW4],
              FW3: [FW4]
              }

# Dependencies={FW1: [FW2],
#               FW2: [FW3],
#               FW3: [FW4]
#               }

WF = Workflow([FW1, FW2, FW3, FW4], Dependencies)

lpad = LaunchPad.auto_load() # loads this based on the FireWorks configuration
lpad.add_wf(WF)

rapidfire(lpad)