#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""
from fireworks import Workflow, LaunchPad
from fireworks.core.rocket_launcher import rapidfire
from triboflow.workflows.subworkflows import MakeAndRelaxSlab_SWF

db_file='/home/mwo/FireWorks/config/db.json'
mp_id = 'mp-81'
functional = 'PBE'
miller = [1, 1, 1]
min_thickness = 10
min_vacuum = 20

WF = MakeAndRelaxSlab_SWF(mp_id, miller, functional)

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
#rapidfire(lpad)