#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""
from fireworks import LaunchPad

from triboflow.workflows.subworkflows import calc_ppes_swf

db_file = "/home/mwo/FireWorks/config/db.json"

WF = calc_ppes_swf(
    interface_name="C001_Ni111_mp-23_mp-48",
    functional="PBE",
    distance_list=[-0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 3.0, 10.0],
    out_name="test_PPES",
    structure_name="unrelaxed_structure",
)

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
# rapidfire(lpad)
