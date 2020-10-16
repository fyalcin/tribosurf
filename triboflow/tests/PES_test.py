#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""
from fireworks import LaunchPad
from fireworks.core.rocket_launcher import rapidfire
from triboflow.workflows.subworkflows import CalcPES_SWF



db_file = '/fs/home/wolloch/git_test/config/db.json'


WF = CalcPES_SWF(interface_name='C001_Ni111_mp-23_mp-48',
                  functional='PBE')


lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
#rapidfire(lpad)
