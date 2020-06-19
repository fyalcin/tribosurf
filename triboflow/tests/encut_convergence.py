#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""
from uuid import uuid4
from fireworks import LaunchPad, Firework, Workflow
from fireworks.core.rocket_launcher import rapidfire
from triboflow.helper_functions import GetLowEnergyStructure, AddBulkToDB
from triboflow.firetasks.utils import FT_PrintFromBulkDB
from triboflow.firetasks.encut_convergence import FT_StartEncutConvo

struct , mp_id = GetLowEnergyStructure('Au')
functional = 'PBE'
db_file='/home/mwo/FireWorks/config/db.json'
AddBulkToDB(struct, mp_id, db_file, functional)

tag = "BM group: {}".format(str(uuid4()))

FT1 = FT_StartEncutConvo(mp_id = mp_id, functional = functional)
FT2 = FT_PrintFromBulkDB(mp_id = mp_id, functional = functional)

FW1 = Firework([FT1], name='Converge Encut Firework')
FW2 = Firework([FT2], name='Print Results Firework')

WF = Workflow([FW1, FW2], {FW1: [FW2]}, name='Encut Convergence Workflow')

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
rapidfire(lpad)
