#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:15:02 2020

@author: mwo
"""
from datetime import datetime

from fireworks import LaunchPad, Workflow, Firework


from fireworks import LaunchPad
from fireworks.core.fworker import FWorker
from fireworks.core.rocket_launcher import rapidfire

from triboflow.workflows.main import heterogeneous_wf_with_surfgen

inputs = {
    "material_1": {"formula": "Si", "mpid": "mp-149"},
    "material_2": {"formula": "Si", "mpid": "mp-149"},
    "computational_params": {
        "functional": "PBE",
        "volume_tolerance": 0.01,
        "BM_tolerance": 0.01,
        "use_vdw": False,
        "use_spin": True,
    },
    "sg_params_1": {
        "slab_thick": 4,
        "vac_thick": 20.0,
        "min_thick_A": 2.0,
        "miller": (1, 0, 0),
    },
    "sg_params_2": {
        "slab_thick": 4,
        "vac_thick": 20.0,
        "min_thick_A": 2.0,
        "miller": (1, 0, 0),
    },
    "sg_filter_1": {"method": "all"},
    "sg_filter_2": {"method": "all"},
    "interface_params": {
        "max_area": 100,
        "max_area_ratio_tol": 0.1,
        "max_length_tol": 0.05,
    },
    "db_file": "/home/yalcin/config_local/db.json",
    "high_level": "surfflow_test",
}

WF = heterogeneous_wf_with_surfgen(inputs)

lpad = LaunchPad.from_file("/home/yalcin/config_local/my_launchpad.yaml")
today = datetime.today().strftime("%Y-%m-%d")
lpad.reset(today)
lpad.add_wf(WF)

# fworker = FWorker.from_file("/home/yalcin/config_local/my_fworker.yaml")
# rapidfire(launchpad=lpad, fworker=fworker, m_dir="/home/yalcin/scratch")
