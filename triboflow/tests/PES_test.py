#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""
from fireworks import LaunchPad
from fireworks.core.rocket_launcher import rapidfire
from triboflow.utils.database import GetDBJSON
from triboflow.workflows.subworkflows import CalcPES_SWF


WF = CalcPES_SWF(interface_name='FeRh001_MgO001_mp-1265_mp-1918',
                  functional='PBE')


lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
#rapidfire(lpad)
