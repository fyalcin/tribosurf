#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""

from fireworks import LaunchPad
from fireworks.core.rocket_launcher import rapidfire
from triboflow.workflows.main import homogeneous_wf


inputs = {
    "material": {
        "formula": "Al",
        "miller": "111",
        "mpid": "mp-134",
        "thick_min": 3,
        "thick_max": 12,
        "thick_incr": 1,
    },
    "computational_params": {
        "functional": "PBE",
        "volume_tolerance": 0.01,
        "BM_tolerance": 0.01,
        "use_vdw": "No",
        "use_spin": "Yes",
        "surfene_thr": 0.01,
        "vacuum": 20,
    },
    "interface_params": {
        "external_pressure": 1.0
    },
}

WF = homogeneous_wf(inputs)

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
# rapidfire(lpad)
