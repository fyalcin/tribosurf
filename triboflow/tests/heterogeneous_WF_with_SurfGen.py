#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""

from fireworks import LaunchPad

from triboflow.workflows.main import heterogeneous_wf

inputs = {
    "material_1": {
        "formula": "Al",
        "mpid": "mp-134"
    },
    "material_2": {
        "formula": "Cu",
        "mpid": "mp-30"
    },
    "computational_params": {
        "functional": "PBE",
        "volume_tolerance": 0.01,
        "BM_tolerance": 0.01,
        "use_vdw": "dftd2",
        "use_spin": True
    },
    "sg_params": {
        "slab_thick": 8,
        "vac_thick": 30.0,
        "min_thick_A": 10.0,
    },
    "interface_params": {
        "max_area": 100,
        "max_area_ratio_tol": 0.1,
        "max_length_tol": 0.05,
    },
}

WF = heterogeneous_wf(inputs)

lpad = LaunchPad.auto_load()
lpad.add_wf(WF)
# rapidfire(lpad)
