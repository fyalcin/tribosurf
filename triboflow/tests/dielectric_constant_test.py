#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 16:29:36 2021

@author: mwo
"""

from fireworks import LaunchPad, Workflow, Firework

from triboflow.firetasks.dielectric import FT_StartDielectric


mpid = 'mp-1265'
functional = 'PBE'

FT = FT_StartDielectric(mp_id=mpid, functional=functional)
WF = Workflow.from_Firework(Firework([FT], name='dielectric_test'))

lpad = LaunchPad().auto_load()
lpad.add_wf(WF)