#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""

from fireworks import LaunchPad

from triboflow.workflows.main import heterogeneous_wf_with_surfgen

inputs = {
    "material_1": {"formula": "Si", "mpid": "mp-149"},
    "material_2": {"formula": "Si", "mpid": "mp-149"},
    "computational_params": {
        "functional": "PBE",
        "volume_tolerance": 0.01,
        "BM_tolerance": 0.01,
        "use_vdw": "dftd2",
        "use_spin": True,
    },
    "sg_params_1": {
        "slab_thick": 6,
        "vac_thick": 20.0,
        "min_thick_A": 4.0,
        "miller": (1, 0, 0)
    },
    "sg_params_2": {
        "slab_thick": 6,
        "vac_thick": 20.0,
        "min_thick_A": 4.0,
        "miller": (1, 0, 0)
    },
    "sg_filter_1": {"method": "all"},
    "sg_filter_2": {"method": "all"},
    "interface_params": {
        "max_area": 100,
        "max_area_ratio_tol": 0.1,
        "max_length_tol": 0.05,
    },
}

WF = heterogeneous_wf_with_surfgen(inputs)

# lpad = LaunchPad.auto_load()
# lpad.add_wf(WF)
# rapidfire(lpad)
