#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 10:37:34 2020

@author: mwo
"""

import os
from fireworks import LaunchPad, Firework, Workflow
from atomate.vasp.fireworks import OptimizeFW
from CommonFiretasks import *
from CommonFireworks import *

indir = os.getcwd()

# create the FireWork
fw1 = FW_ReadWFInputs('Input_test.in')

fw2 = Firework(FT_FetchStructureFromInput(input_dict_name='Input_test.in'))

fw3 = Firework(FT_SpawnOptimizeFW())

fw4 = Firework(FT_PrintSpec())

wf = Workflow([fw1, fw2, fw3, fw4], {fw1: [fw2], fw2: [fw3], fw3: [fw4]})

# finally, instatiate the LaunchPad and add the workflow to it
lpad = LaunchPad.auto_load() # loads this based on the FireWorks configuration
lpad.add_wf(wf)
